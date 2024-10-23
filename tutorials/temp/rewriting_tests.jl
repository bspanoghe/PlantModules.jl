using BenchmarkTools, Infiltrator, Revise
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots

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

rule = Rule(
    Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchingrule
)

axiom = Stem([0.5, 5.0]) + (Leaf([3.0, 1.0, 0.1]), Leaf([5.0, 3.0, 0.1]))

plant = Graph(axiom = axiom, rules = (rule,))
num_iterations = 3
for _ in 1:num_iterations
    rewrite!(plant)
end

convert_to_MTG(plant) |> PlantModules.MultiScaleTreeGraph.DataFrame

## Environment

soil_graph = Soil()
air_graph = Air()

## Combined

graphs = [plant, soil_graph, air_graph]
intergraph_connections = [(1, 2) => (PlantModules.root(plant), :Soil), (1, 3) => (:Leaf, :Air)]
struct_connections = PlantStructure(graphs, intergraph_connections)

# Functional processes

ϵ = 0.03
ϕ = 3e-3
module_defaults = Dict(
	:Stem => Dict(:shape => Cilinder(ϵ_D = fill(ϵ, 2), ϕ_D = fill(ϕ, 2)),
        :M => 400e-6),
	:Leaf => Dict(:shape => Cuboid(ϵ_D = fill(ϵ, 3), ϕ_D = fill(ϕ, 3)),
        :M => 450e-6, :K_s => 1e-4),
	:Soil => Dict(:W_max => 1e3, :T => 288.15),
	:Air => Dict(:W_r => 0.8)
)

connecting_modules = [
	(:Soil, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Leaf) => (const_hydraulic_connection, Dict()),
	(:Leaf, :Air) => (evaporation_connection, Dict()),
]

func_connections = PlantFunctionality(; module_defaults, connecting_modules)

# Coupling

module_coupling = Dict(
    :Stem => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Leaf => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Soil => [environmental_module, Ψ_soil_module, constant_K_module],
    :Air => [environmental_module, Ψ_air_module],
)

# Rev her up

system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false);
sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24));
# @btime sol = solve(prob);
sol = solve(prob);

#=
SOLTIMES

3 (45): 0.1s
4: 1.8s
5 (189): 6s
6 (381): 28s
7 (765): 135s

PARAM INFLUENCES (3 rewrites)
Base (ϵ = 3, ϕ = 3e-3, K = 1000 (Leaf: 1e-4), W_max soil: 1e5): 104ms

elastic modulus ϵ (for very low values => no growth; same results for medium small to large values)
    - 0.0003: 51ms
    - 0.003: 44ms
    - 0.03: 49ms
    - 0.3: 83ms
    ------------
    - 30: 102ms
    - 300: 116ms
    - 3000: 110ms
    - 300_000: 155ms

extensibility ϕ (drastically different simulation results)
    - 3e-5: 104ms
    - 3e-4: 95ms
    -----------
    - 3e-2: 166ms
    - 3e-1: 196ms (soilwater halved)

evaporation change:
- smooth: 110-120ms
- abrupt: 150-170ms
=#

histogram(sol.t)

# Plotting

plotgraph(sol, graphs[1], varname = :W)
plotgraph(sol, graphs[2], varname = :W)

plotnode(sol, PlantModules.root(graphs[1]), varname = :D)
plotnode(sol, PlantModules.root(graphs[1]), varname = :D, ylims = (0, 1))
plotnode(sol, PlantModules.nodes(graphs[1])[end], varname = :D)

plotgraph(sol, graphs[1:2], varname = :Ψ)
plotgraph(sol, graphs[1], varname = :D)
plotgraph(sol, graphs[1], varname = :P, structmod = :Stem)
plotgraph(sol, graphs[1], varname = :Π)

plotgraph(sol, graphs[3], varname = :ΣF)