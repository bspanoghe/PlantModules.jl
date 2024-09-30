using Revise, BenchmarkTools, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, MultiScaleTreeGraph
using ModelingToolkit, OrdinaryDiffEq, Unitful
using PlantBiophysics, PlantBiophysics.PlantMeteo#, PlantSimEngine
using Memoization
using Plots; import GLMakie.draw

# Structural modules #

## Plant
plant_graph = readXEG("tutorials/temp/structures/beech10.xeg") #! change to beech
# convert_to_PG(plant_graph) |> draw

function get_mtg()
	plant_graph = readXEG("tutorials/temp/structures/beech.xeg")
	mtg = convert_to_MTG(plant_graph)

	### Setting attributes right
	DataFrame(mtg, [:diameter, :length, :width])
	
	function combine_dimensions(l, d, w)
		if all(isnothing.([l, d, w]))
			return nothing
		elseif isnothing(w)
			return 1e2*[d/2, l]
		else
			return 1e2*[l, w, 5e-4]
		end
	end
	
	transform!(mtg, [:length, :diameter, :width] => combine_dimensions => :D)
	DataFrame(mtg, [:D])
	
	### Inspecting what kind of structural modules are in here
	# me_structmods = [PlantModules.structmod(node) for node in PlantModules.nodes(mtg)] |> unique
	
	# for me_structmod in me_structmods
	# 	dimensions = [node_attributes(node)[:D] for node in PlantModules.nodes(mtg) if PlantModules.structmod(node) == me_structmod]
	# 	num_nodes = length(dimensions)
	# 	println("There are $num_nodes nodes of type $me_structmod $( all(ismissing.(dimensions)) ? "(dimensions undefined)" : "")")
	# end
	
	traverse!(mtg, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")
	mtg = delete_nodes!(mtg, filter_fun = node -> isnothing(node.D))

	# descendants(mtg, symbol = "Internode", self = true) |> DataFrame
	# descendants(mtg, symbol = "Shoot", self = true) |> DataFrame
	# descendants(mtg, symbol = "Leaf", self = true) |> DataFrame
end

# mtg = delete_nodes!(mtg, symbol = "Shoot")

function prunetree!(mtg, target_length; remaining_length = length(mtg))
	if remaining_length - length(mtg) > target_length
		branchnodes = descendants(mtg, self = true)
		delete_nodes!(mtg, filter_fun = node -> node in branchnodes)
	end

	for chnode in children(mtg)
		prunetree!(chnode, target_length, remaining_length = length(mtg))
	end
end

function create_MWE()
	mtg = get_mtg()
	[prunetree!(mtg, 60) for _ in 1:2]; # no work
	# mtg = delete_node!(mtg)
	# mtg = delete_node!(mtg)
	# mtg = delete_node!(mtg)

	# lil_branch = descendants(mtg[1][1][1][1][1][1][1][1][1][1][1][1], self = true)
	# mtg = delete_nodes!(mtg, filter_fun = node -> node in lil_branch)

	# bigger_branch = descendants(mtg[2][1][1][1][1][1][1][1][1], self = true)
	# mtg = delete_nodes!(mtg, filter_fun = node -> node in bigger_branch)

	# shoots = descendants(mtg, symbol = "Shoot", self = true)
	# mtg = delete_nodes!(mtg, filter_fun = node -> node in shoots)

	return mtg
end

mtg = create_MWE()
length(mtg)
# length(mtg)
# convert_to_PG(mtg) |> draw

# mtg = get_mtg()
# [prunetree!(mtg, 40) for _ in 1:3]; # yes bueno

# DataFrame(mtg, [:D])
# branchingpoint = mtg[1][1][1]
# DataFrame(branchingpoint, [:D])
# DataFrame(branchingpoint[1], [:D])
# DataFrame(branchingpoint[2], [:D])

# supervolume(mtg) = volume(length(mtg.D) == 2 ? Cilinder(ϵ_D = [0, 0], ϕ_D = [0, 0]) : Cuboid(ϵ_D = [0, 0, 0], ϕ_D = [0, 0, 0]), mtg.D)

# minimum(supervolume.(descendants(mtg, self = true)))
# mean(supervolume.(descendants(mtg, self = true)))


# function checkthemconnections(mtg)
# 	vol = supervolume(mtg)
# 	chvols = [supervolume(chnode) for chnode in PlantModules.children(mtg, mtg)]
# 	println("$(round.(chvols ./ vol, digits = 2))")

# 	for chnode in PlantModules.children(mtg, mtg)
# 		checkthemconnections(chnode)
# 	end
# end

# checkthemconnections(mtg)

#! 

# import PlantGraphs: Node
# mutable struct Internode <: Node
# 	D::Vector
# end
# mutable struct Leaf <: Node
# 	D::Vector
# end
# D_internode

# mtg = InterNode()

## Environment

struct Soil <: PlantGraphs.Node end
struct Air <: PlantGraphs.Node end

soil_graph = Soil()
air_graph = Air()

graphs = [mtg, soil_graph, air_graph]

## connections

intergraph_connections = [[1, 2] => (mtg, :Soil), [1, 3] => (:Leaf, :Air)]
struct_connections = PlantStructure(graphs, intergraph_connections)

# Functional modules #

## New functional modules

# PAR_data = [max(0, 400 * sin(t/24*2*pi - 8)) + randn() for t in 0:10*24]
# @memoize get_PAR_flux(t) = PAR_data[floor(Int, t+1)] + (t-floor(t)) * PAR_data[ceil(Int, t+1)]
get_PAR_flux(t) = max(0, 400 * sin(t/24*2*pi - 8))
@register_symbolic get_PAR_flux(t)

@memoize function get_assimilation_rate(PAR_flux, T, LAI, k)
	Kelvin_to_C = -273.15
	meteo = Atmosphere(T = T + Kelvin_to_C, Wind = 1.0, P = 101.3, Rh = 0.65, Ri_PAR_f = PAR_flux)
	m = ModelList(
		Fvcb(), # calculate CO2 assimilation rate
		Medlyn(0.03, 0.92), # calculate stomatal conductance, see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2486.2010.02375.x
		Beer(k), # calculate amount of light intercepted
		status = (Tₗ = meteo[:T], LAI = LAI, Cₛ = meteo[:Cₐ], Dₗ = meteo[:VPD], RI_PAR_f = meteo[:Ri_PAR_f])
	)
	run!(m, meteo)
	return only(m[:A]) |> x -> max(x, 0) # extract result of the first (and only) timestep
end

@register_symbolic get_assimilation_rate(PAR_flux, T, LAI, k)

import .PlantModules: t, d

function photosynthesis_module(; name, T, M, shape)
	@constants (
		uc1 = (10^-6 * 10^-4 * 60^2), [description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)", unit = u"(mol/cm^2/hr) / (µmol/m^2/s)"],
	)
	@parameters (
		T = T, [description = "Temperature", unit = u"K"],
		LAI = 2.0, [description = "Leaf Area Index", unit = u"cm^2 / cm^2"],
		k = 0.5, [description = "Light extinction coefficient", unit = u"N/N"],
		carbon_decay_rate = 0.1, [description = "Rate at which carbon is consumed for growth", unit = u"hr^-1"],
	)

	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
		PF(t), [description = "Incoming PAR flux", unit = u"J / s / m^2"], #! make sure not to use variable name from other func mod used in same struct mod
		A(t), [description = "Carbon assimilation rate", unit = u"µmol / m^2 / s"],
		D(t)[1:length(shape.ϵ_D)], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		PF ~ get_PAR_flux(t)
		A ~ get_assimilation_rate(PF, T, LAI, k)
        d(M) ~ uc1 * A * cross_area(shape, D) / volume(shape, D) - carbon_decay_rate*M # convert µmol => mol and s^-1 => hr^-1
		#! change carbon decay rate into maintenance term (~ compartment size) and growth term (~ compartment growth)
		#! add buffer term? (#starch)
    ]
    return ODESystem(eqs, t; name, checks = false) #! checks back to true?
end

function waterdependent_K_module(; name, K_max)
	@parameters (
		K_max = K_max, [description = "Hydraulic conductivity of compartment at maximum water content", unit = u"g / hr / MPa"],
	)
    @variables (
		W_r(t), [description = "Relative water content of compartment", unit = u"g / g"],
		K(t), [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
    )

    eqs = [
		K ~ K_max * exp(-5*(1-W_r)^3)
    ]

    return ODESystem(eqs, t; name)
end

C_stem = 300e-6
C_shoot = 350e-6
C_leaf = 400e-6

extra_defaults = Dict( 
	photosynthesis_module => Dict(:shape => Cuboid(ϵ_D = [0, 0, 0], ϕ_D = [0, 0, 0]), :T => 293.15, :M => 0),
	waterdependent_K_module => Dict(:K_max => 10), 
)

module_defaults = Dict(
	:Internode => Dict(:shape => Cilinder(ϵ_D = [10.0, 300.0], ϕ_D = [1e-3, 1e-5]), :M => C_stem),
	:Shoot => Dict(:shape => Cilinder(ϵ_D = [10.0, 300.0], ϕ_D = [3e-3, 3e-5]), :M => C_shoot),
	:Leaf => Dict(:shape => Cuboid(ϵ_D = [15.0, 10.0, 1000.0], ϕ_D = [3e-3, 3e-3, 1e-5]), :M => C_leaf, :K_s => 0.01),
	:Soil => Dict(:W_max => 10000.0, :T => 293.15, :W_r => 0.9),
	:Air => Dict(:W_r => 0.8)
)

connecting_modules = [
	(:Soil, :Internode) => (hydraulic_connection, Dict()),
    (:Internode, :Internode) => (hydraulic_connection, Dict()),
	(:Internode, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Leaf) => (const_hydraulic_connection, Dict(:K => 10)),
	(:Internode, :Leaf) => (hydraulic_connection, Dict()),
    (:Leaf, :Air) => (hydraulic_connection, Dict())
]

func_connections = PlantFunctionality(module_defaults = module_defaults, 
	connecting_modules = connecting_modules, extra_defaults = extra_defaults
)

## Connect them to structure

module_coupling = Dict( #! photosynthesis_module for leaf
	:Internode => [hydraulic_module, constant_carbon_module, sizedep_K_module],
	:Shoot => [hydraulic_module, constant_carbon_module, sizedep_K_module],
	:Leaf => [hydraulic_module, constant_carbon_module, sizedep_K_module],
	:Soil => [environmental_module, Ψ_soil_module, waterdependent_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
)

# Gettem #

system = PlantModules.generate_system(struct_connections, func_connections, module_coupling, checkunits = false)

sys_simpl = structural_simplify(system)
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
@time sol = solve(prob);
# @btime sol = solve(prob);

PlantModules.plotgraph(sol, graphs[1], func_varname = :W, struct_module = :Leaf)
PlantModules.plotgraph(sol, graphs[1], func_varname = :W, struct_module = :Internode)

PlantModules.plotgraph(sol, graphs[1], func_varname = :Ψ, struct_module = :Internode)


PlantModules.plotgraph(sol, graphs[2], func_varname = :W)
PlantModules.plotgraph(sol, graphs[3], func_varname = :W)

PlantModules.plotgraph(sol, graphs[1], func_varname = :ΣF)
PlantModules.plotgraph(sol, graphs[2], func_varname = :ΣF)
PlantModules.plotgraph(sol, graphs[3], func_varname = :ΣF)

PlantModules.plotgraph(sol, graphs[1], func_varname = :M)
PlantModules.plotgraph(sol, graphs[1:2], func_varname = :Ψ)
PlantModules.plotgraph(sol, graphs[1], func_varname = :A)
PlantModules.plotgraph(sol, graphs[1], func_varname = :PF)

PlantModules.plotgraph(sol, graphs[1], func_varname = :P)
PlantModules.plotgraph(sol, graphs[1], func_varname = :D, struct_module = :Leaf)


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