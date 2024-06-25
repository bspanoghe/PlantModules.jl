using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs, MultiScaleTreeGraph
import PlantGraphs: Node
using ModelingToolkit, DifferentialEquations, Unitful
using Plots; import GLMakie.draw
using BenchmarkTools

mutable struct Root <: Node end
mutable struct Stem <: Node end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

C_root = 300e-6 # (mol/cm^3 = 10^6 µmol/g) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 400e-6 # Idem here
C_leaf = 450e-6 # And here!
# loosely based on https://www.researchgate.net/figure/Total-soluble-sugar-content-in-the-leaf-a-and-root-tissues-b-of-Homjan-HJ_fig3_226164026

PlantModules.default_params
PlantModules.default_u0s

default_changes = (Γ = 0.4, P = 0.2, T = 293.15)

default_params, default_u0s = PlantModules.alter_defaults(default_changes)

module_defaults = (
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 10.0, 20.0], ϕ_D = [0.7, 0.1, 0.05]), D = [30, 3, 1], M = C_root), #! changed C to M_const #! include check whether module_defaults variabkle names correspond with default_params / default_u0s ?
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 15.0], ϕ_D = [0.8, 0.03]), D = [1.5, 10], M = C_stem),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 500.0, T = 288.15),
	Air = (W_r = 0.8,)
)

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

if !isdefined(Main, :plant_graph) # for debugging mostly
	plant_graph = Root() + Stem() + (Leaf([2]), Leaf([2.5]))
end

# draw(plant_graph, resolution = (500, 400))

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = [graphs, intergraph_connections]

connecting_modules = [
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 8]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 4]), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Stem, :Leaf) => (PlantModules.hydraulic_connection, [:K => 4]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-3]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-3]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

func_connections = [connecting_modules, PlantModules.multi_connection_eqs]

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))	 #! Generate warning for missing variable defaults that aren't dummies
@time sol = solve(prob);

PlantModules.plotgraph(sol, graphs[1]) .|> display
PlantModules.plotgraph(sol, graphs[2]) .|> display
PlantModules.plotgraph(sol, graphs[3]) .|> display
PlantModules.plotgraph(sol, graphs[1:2], func_varname = :Ψ)
PlantModules.plotgraph(sol, graphs[1:2], func_varname = :W)
