using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs, MultiScaleTreeGraph
import PlantGraphs: Node
using ModelingToolkit, DifferentialEquations, Unitful
using Plots; import GLMakie.draw
using PlantBiophysics, PlantBiophysics.PlantMeteo, PlantSimEngine
using Memoization

# Structural modules #

mutable struct Root <: Node end
mutable struct Internode <: Node end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

if !isdefined(Main, :plant_graph) # for debugging mostly
	plant_graph = Root() + Internode() + (Internode() + Leaf([5, 3, 0.08]), Internode() + (Leaf([4, 1, 0.05]), Leaf([2, 0.3, 0.04])), Leaf([8, 4, 0.1]))
end

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air)]

struct_connections = [graphs, intergraph_connections]

# Functional modules #

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
leafarea(::Cuboid, D::AbstractArray) = D[1] * D[2]

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
        d(M) ~ uc1 * A * leafarea(shape, D) / volume(shape, D) - carbon_decay_rate*M # convert µmol => mol and s^-1 => hr^-1
    ]
    return ODESystem(eqs, t; name, checks = false) #! checks back to true?
end

## Connections

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Internode, :Leaf],
    photosynthesis_module => [:Leaf],
	PlantModules.constant_carbon_module => [:Root, :Internode],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air],
]

connecting_modules = [
    (:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 8]),
	(:Root, :Internode) => (PlantModules.hydraulic_connection, [:K => 10]),
    (:Internode, :Internode) => (PlantModules.hydraulic_connection, [:K => 2]),
	(:Internode, :Leaf) => (PlantModules.hydraulic_connection, [:K => 1]),
    (:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-3])
]
func_connections = [connecting_modules, PlantModules.multi_connection_eqs]


## Parameters

default_changes = (Γ = 0.4, P = 0.2, D = [0.5, 5.0])
default_params, default_u0s = PlantModules.alter_defaults(default_changes)

default_params = merge(default_params, 
	(photosynthesis_module = (shape = Cuboid(ϵ_D = [0, 0, 0], ϕ_D = [0, 0, 0]), T = 293.15),),
)

default_u0s = merge(default_u0s,
	(photosynthesis_module = (M = 0,),),
)

C_root = 300e-6
C_stem = 400e-6 
C_leaf = 450e-6 

module_defaults = (
    Root = (shape = Cilinder(ϵ_D = [10.0, 300.0], ϕ_D = [3e-3, 3e-5]), M = C_root),
	Internode = (shape = Cilinder(ϵ_D = [10.0, 300.0], ϕ_D = [1e-3, 1e-5]), M = C_stem),
	Leaf = (shape = Cuboid(ϵ_D = [5.0, 3.0, 500.0], ϕ_D = [5e-3, 3e-3, 1e-5]), M = C_leaf),
	Soil = (W_max = 10000.0, T = 293.15, W_r = 0.9),
	Air = (W_r = 0.8,)
)

# GO #

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))	 #! Generate warning for missing variable defaults that aren't dummies
@time sol = solve(prob);


plotgraph(sol, graphs[1], func_varname = :M, struct_module = :Leaf);
plot!(yformatter = y -> y * 1e3);
ylabel!("active metabolites (1e³ mol / cm³)")
savefig("./plots/Leaf_M.svg")

plotgraph(sol, graphs[1], func_varname = :W)
ylabel!("water (g)")
savefig("./plots/Plant_W.svg")
