using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs, MultiScaleTreeGraph
using ModelingToolkit, DifferentialEquations, Unitful
using PlantBiophysics, PlantBiophysics.PlantMeteo, PlantSimEngine
using Surrogates
using Plots; import GLMakie.draw

# Structural modules #

## Plant
plant_graph = readXEG("tutorials/temp/structures/beech0.xeg") #! change to beech
# convert_to_PG(plant_graph) |> draw

mtg = convert_to_MTG(plant_graph)

### Setting attributes right
DataFrame(mtg, [:diameter, :length, :width])

function combine_dimensions(l, d, w)
	if all(isnothing.([l, d, w]))
		return nothing
	elseif isnothing(w)
		return [l, d]
	else
		return [l, w, 1.5e-4]
	end
end

transform!(mtg, [:length, :diameter, :width] => combine_dimensions => :D)
DataFrame(mtg, [:D])

function get_shape(D)
	if isnothing(D)
		return nothing
	elseif length(D) == 2
		return Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.1, 0.005])
	else
		Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.1, 0.01, 0.005])
	end
end
transform!(mtg, :D => get_shape => :shape)
DataFrame(mtg, [:shape])

### Inspecting what kind of structural modules are in here
me_structmods = [PlantModules.structmod(node) for node in PlantModules.nodes(mtg)] |> unique

for me_structmod in me_structmods
	dimensions = [node_attributes(node)[:D] for node in PlantModules.nodes(mtg) if PlantModules.structmod(node) == me_structmod]
	num_nodes = length(dimensions)
	println("There are $num_nodes nodes of type $me_structmod $( all(ismissing.(dimensions)) ? "(dimensions undefined)" : "")")
end

traverse!(mtg, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")
mtg = delete_nodes!(mtg, filter_fun = node -> isnothing(node.D))

mtg |> DataFrame

descendants(mtg, symbol = "Internode", self = true) |> DataFrame
descendants(mtg, symbol = "Shoot", self = true) |> DataFrame
descendants(mtg, symbol = "Leaf", self = true) |> DataFrame

## Environment

struct Soil <: PlantGraphs.Node end
struct Air <: PlantGraphs.Node end

soil_graph = Soil()
air_graph = Air()

graphs = [mtg, soil_graph, air_graph]

## connections

intergraph_connections = [[1, 2] => (mtg, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)]
struct_connections = [graphs, intergraph_connections]

# Functional modules #

## New functional modules

PAR_data = [max(0, 400 * sin(t/24*2*pi - 8)) + randn() for t in 0:10*24]
get_PAR_flux(t) = PAR_data[floor(Int, t+1)] + (t-floor(t)) * PAR_data[ceil(Int, t+1)]
@register_symbolic get_PAR_flux(t)

function get_assimilation_rate(PAR_flux, T, LAI, k)
	Kelvin_to_C = -273.15
	meteo = Atmosphere(T = T + Kelvin_to_C, Wind = 1.0, P = 101.3, Rh = 0.65, Ri_PAR_f = PAR_flux)
	m = ModelList(
		Fvcb(), # calculate CO2 assimilation rate
		Medlyn(0.03, 0.92), # calculate stomatal conductance, see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2486.2010.02375.x
		Beer(k), # calculate amount of light intercepted
		status = (Tₗ = meteo[:T], LAI = LAI, Cₛ = meteo[:Cₐ], Dₗ = meteo[:VPD], RI_PAR_f = meteo[:Ri_PAR_f])
	)
	run!(m, meteo)
	return only(m[:A]) # extract result of the first (and only) timestep
end

n_samples = 1000
lower_bound = [0.0, 273.15]
upper_bound = [600.0, 303.15]
PAR_Ts = sample(n_samples, lower_bound, upper_bound, RandomSample())
As = get_assimilation_rate.(first.(PAR_Ts), last.(PAR_Ts), 2.0, 0.5)
surr = SecondOrderPolynomialSurrogate(PAR_Ts, As, lower_bound, upper_bound)

test_PAR_Ts = sample(n_samples, lower_bound, upper_bound, RandomSample())
test_As = get_assimilation_rate.(first.(test_PAR_Ts), last.(test_PAR_Ts), 2.0, 0.5)
mean((test_As - surr.(test_PAR_Ts)).^2)
@btime surr([400.2, 295.3])

get_est_assimilation_rate(PAR, T) = [surr([PAR_i, T]) for PAR_i in PAR]
@register_symbolic get_est_assimilation_rate(PAR, T)

leafarea(::Cuboid, D::AbstractArray) = D[1] * D[2]

import .PlantModules: t, d

function photosynthesis_module(; name, T, M, shape, D)
	@constants (
		uc1 = (10^-6 * 60^2), [description = "Unit conversion from (µmol/s) to (mol/hr)", unit = u"(mol/hr) / (µmol/s)"],
	)
	@parameters (
		T = T, [description = "Temperature", unit = u"K"],
		LAI = 2.0, [description = "Leaf Area Index", unit = u"m^2 / m^2"],
		k = 0.5, [description = "Light extinction coefficient", unit = u"N/N"],
		carbon_decay_rate = 0.1, [description = "Rate at which carbon is consumed for growth", unit = u"hr^-1"],
	)

	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / m^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
		PF(t), [description = "Incoming PAR flux", unit = u"J / s / m^2"], #! make sure not to use variable name from other func mod used in same struct mod
		A(t), [description = "Carbon assimilation rate", unit = u"µmol / m^2 / s"],
    )

    eqs = [
		PF ~ get_PAR_flux(t)
		A ~ get_est_assimilation_rate(PAR, T)
        d(M) ~ uc1 * leafarea(shape, D) * A / volume(shape, D) - carbon_decay_rate*M # convert µmol => mol and s^-1 => hr^-1
    ]
    return ODESystem(eqs, t; name, checks = false) #! checks back to true?
end

## Connect them to structure

module_coupling = [
	photosynthesis_module => [:Leaf],
	PlantModules.hydraulic_module => [:Internode, :Shoot, :Leaf],
	PlantModules.constant_carbon_module => [:Internode, :Shoot],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air],
]

connecting_modules = [
	() => PlantModules.hydraulic_connection,
	(:Soil, :Internode) => (PlantModules.hydraulic_connection, [:K => 3]),
    (:Internode, :Internode) => (PlantModules.hydraulic_connection, [:K => 3]),
	(:Internode, :Shoot) => (PlantModules.hydraulic_connection, [:K => 2]),
	(:Internode, :Leaf) => (PlantModules.hydraulic_connection, [:K => 2]),
    (:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 0.01])
]

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

func_connections = [connecting_modules, get_connection_eqs]

## Tweak parameters

C_stem = 300
C_shoot = 350
C_leaf = 400

default_params = merge(PlantModules.default_params, 
	(photosynthesis_module = (
		T = 298.15,	
		shape = Cuboid(ϵ_D = [0, 0, 0], ϕ_D = [0, 0, 0]),
	),)
)

default_u0s = merge(PlantModules.default_u0s,
	(photosynthesis_module = (
		M = C_leaf,
		D = [0, 0, 0]
	),)
)

module_defaults = (
	Internode = (shape = Cilinder(ϵ_D = [5.0, 0.3], ϕ_D = [0.1, 0.01]), M = C_stem),
	Shoot = (shape = Cilinder(ϵ_D = [5.0, 0.3], ϕ_D = [0.1, 0.01]), M = C_shoot),
	Leaf = (shape = Cuboid(ϵ_D = [0.5, 0.5, 0.01], ϕ_D = [0.1, 0.1, 0.01]), M = C_leaf),
	Soil = (W_max = 100000.0, T = 293.15),
	Air = ()
)

# Gettem #

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system)
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
sol = solve(prob)

PlantModules.plotgraph(sol, graphs[1], func_varname = :W)
PlantModules.plotgraph(sol, graphs[2], func_varname = :W)
PlantModules.plotgraph(sol, graphs[1], func_varname = :M)
PlantModules.plotgraph(sol, graphs[1], func_varname = :Ψ)
PlantModules.plotgraph(sol, graphs[1:2], func_varname = :Ψ)
PlantModules.plotgraph(sol, graphs[1], func_varname = :A)

# plotspeed #

xs = [rand(100) for _ in 1:1000]
ys = [rand(100) for _ in 1:1000]

function plot1(xs, ys)
    plt = plot()
    for (x, y) in zip(xs, ys)
        plot!(plt, x, y)
    end
    plt
end

function plot2(xs, ys)
    plot(reduce(vcat, xs), reduce(vcat, ys))
end

@btime plot1($xs, $ys) ; # 431.335 ms (2231129 allocations: 121.87 MiB)
@btime plot2($xs, $ys) ; # 2.844 ms (2530 allocations: 3.20 MiB)

xs = [i%2 == 0 ? NaN : 1:10 for i in 1:6]
ys = [i%2 == 0 ? NaN : rand(10) for i in 1:6]
plot2(xs, ys)