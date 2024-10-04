using BenchmarkTools, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots; import GLMakie.draw


# Structure

## Plant

mutable struct Stem <: Node
	D::Vector
end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

function branchingrule(x)
    leaf_D = [2.0, 1.0, 0.1]
    parent_D = data(x).D
    child_D = [parent_D[1] / sqrt(2), 1.0]
    return Stem(parent_D) + (Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)), Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)))
end

rule = Rule(Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchingrule
)

axiom = Stem([0.5, 5.0]) + (Leaf([3.0, 1.0, 0.1]), Leaf([5.0, 3.0, 0.1]))

plant = Graph(axiom = axiom, rules = (rule,))
num_iterations = 4
for _ in 1:num_iterations
    rewrite!(plant)
end

convert_to_MTG(plant) |> PlantModules.MultiScaleTreeGraph.DataFrame

## Environment

soil_graph = Soil()
air_graph = Air()

## Combined

graphs = [plant, soil_graph, air_graph]
intergraph_connections = [(1, 2) => (PlantModules.root(plant), :Soil), (1, 3) => (:Leaf, :Air), (2, 3) => (:Soil, :Air)]
struct_connections = PlantStructure(graphs, intergraph_connections)


# Functional processes

module_defaults = Dict(
	:Stem => Dict(
        :shape => Cilinder(ϵ_D = [2.0, 4.5], ϕ_D = 1e-3 .* [8, 3]), :D => [1.5, 10], :M => 400e-6, :K_s => 10,
        :P => -0.3 + 8.314*298.15*400e-6
        ),
	:Leaf => Dict(
        :shape => Cuboid(ϵ_D = [1.5, 1.5, 10.0], ϕ_D = 1e-3 .* [3, 3, 0.1]), :M => 450e-6, :K_s => 1e-5,
        :P => -0.3 + 8.314*298.15*450e-6
        ),
	:Soil => Dict(:W_max => 10000.0, :T => 288.15, :K => 10),
	:Air => Dict(:W_r => 0.8, :K => 0.1)
)

connecting_modules = [
	(:Soil, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Leaf) => (const_hydraulic_connection, Dict(:K => 10)),
	(:Leaf, :Air) => (hydraulic_connection, Dict()),
	(:Soil, :Air) => (hydraulic_connection, Dict())
]

func_connections = PlantFunctionality(module_defaults = module_defaults, connecting_modules = connecting_modules)

# Coupling 

module_coupling = Dict(
    :Stem => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Leaf => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Soil => [environmental_module, Ψ_soil_module, constant_K_module],
    :Air => [environmental_module, Ψ_air_module, constant_K_module],
)

# Rev her up

system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false)
sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
@time sol = solve(prob);

# Plotting

plotgraph(sol, graphs[1], func_varname = :W)
plotgraph(sol, graphs[1], func_varname = :W, ylims = (0, 1))
plotgraph(sol, graphs[2], func_varname = :W)

plotgraph(sol, graphs[1], func_varname = :Ψ, struct_module = :Stem)
plotgraph(sol, graphs[1], func_varname = :P, struct_module = :Stem)
plotgraph(sol, graphs[1], func_varname = :Π, struct_module = :Stem)

plotgraph(sol, graphs[2], func_varname = :Ψ)