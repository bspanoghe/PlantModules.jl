using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs
using ModelingToolkit, DifferentialEquations, Unitful
using Plots; import GLMakie.draw
using BenchmarkTools

mutable struct Stem <: Node
	D::Vector
end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

axiom = Stem([0.5, 5.0]) + (Leaf([3.0, 1.0, 0.1]), Leaf([5.0, 3.0, 0.1]))
function branchingrule(x)
    leaf_D = [3.0, 1.0, 0.1]
    parent_D = data(x).D
    child_D = [parent_D[1] / sqrt(2), 1.0]
    return Stem(parent_D) + (Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)), Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)))
end

rule1 = Rule(Stem, rhs = branchingrule)
plant = Graph(axiom = axiom, rules = (rule1,))

for _ in 1:3
    rewrite!(plant)
end

draw(plant)


C_root = 300e-6
C_stem = 400e-6
C_leaf = 450e-6

module_defaults = (
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