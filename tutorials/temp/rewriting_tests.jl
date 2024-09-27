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

function branchlessrule(x)
    parent_D = data(x).D
    return Stem(parent_D) + Stem(parent_D)
end

function growify!(plant, growthrate)
    for node in apply(plant, Query(Stem))
        node.D = [sqrt(growthrate)*node.D[1], growthrate*node.D[2]]
    end

	for node in apply(plant, Query(Leaf))
        node.D = [growthrate*node.D[1], growthrate*node.D[2], node.D[3]]
    end
end

rule1 = Rule(Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchingrule
)

rule2 = Rule(Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchlessrule
)

axiom = Stem([0.5, 5.0]) + (Leaf([3.0, 1.0, 0.1]), Leaf([5.0, 3.0, 0.1]))

plant = Graph(axiom = axiom, rules = (rule1,))
num_iterations = 7
for _ in 1:num_iterations
    rewrite!(plant)
	growify!(plant, 1.03)
end

# draw(plant)
convert_to_MTG(plant) |> PlantModules.MultiScaleTreeGraph.DataFrame

## Environment

soil_graph = Soil()
air_graph = Air()

## Combined

graphs = [plant, soil_graph, air_graph]
intergraph_connections = [(1, 2) => (PlantModules.root(plant), :Soil), (1, 3) => (:Leaf, :Air), (2, 3) => (:Soil, :Air)]
struct_connections = PlantStructure(graphs, intergraph_connections)


# Functional processes

t = PlantModules.t

function K_surface_module(; name, K_surface, shape::PlantModules.Shape)
    num_D = length(shape.ϵ_D)

    @parameters (
        K_s = K_surface, [description = "Specific hydraulic conductivity of compartment", unit = u"g / hr / MPa / cm^2"],
    )
    @variables (
        K_surface(t), [description = "Hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
		D(t)[1:num_D], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		K_surface ~ K_s * surface_area(shape, D)/2
    ]

    return ODESystem(eqs, t; name)
end

function surface_hydraulic_connection(; name)
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1", unit = u"g / hr / MPa"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2", unit = u"g / hr / MPa"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ min(K_1, K_2) * (Ψ_2 - Ψ_1) #! mean?
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order) = 
        if !reverse_order
            [
                connection_MTK.Ψ_1 ~ node_MTK.Ψ,
                connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
                connection_MTK.K_1 ~ node_MTK.K_surface,
                connection_MTK.K_2 ~ nb_node_MTK.K,
            ]
        else
            [
                connection_MTK.Ψ_1 ~ node_MTK.Ψ,
                connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
                connection_MTK.K_1 ~ node_MTK.K,
                connection_MTK.K_2 ~ nb_node_MTK.K_surface,
            ]
        end

    return ODESystem(eqs, t; name), get_connection_eqset
end

C_root = 300e-6
C_stem = 400e-6
C_leaf = 450e-6

extra_defaults = Dict(
    K_surface_module => Dict(:K_surface => 0, :shape => Sphere())
)

module_defaults = Dict(
	:Stem => Dict(:shape => Cilinder(ϵ_D = [2.0, 4.5], ϕ_D = 1e-3 .* [8, 3]), :D => [1.5, 10], :M => C_stem, :K_s => 10),
	:Leaf => Dict(:shape => Cuboid(ϵ_D = [1.5, 1.5, 10.0], ϕ_D = 1e-3 .* [3, 3, 0.1]), :M => C_leaf, :K => 10, :K_surface => 1e-5),
	:Soil => Dict(:W_max => 10000.0, :T => 288.15, :K => 10),
	:Air => Dict(:W_r => 0.8, :K => 0.1)
)

connecting_modules = [
	(:Soil, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Leaf) => (hydraulic_connection, Dict()),
	(:Leaf, :Air) => (surface_hydraulic_connection, Dict()),
	(:Soil, :Air) => (hydraulic_connection, Dict())
]

func_connections = PlantFunctionality(module_defaults = module_defaults, connecting_modules = connecting_modules, extra_defaults = extra_defaults)

# Coupling 

module_coupling = Dict(
    :Stem => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Leaf => [hydraulic_module, constant_carbon_module, constant_K_module, K_surface_module],
    :Soil => [environmental_module, Ψ_soil_module, constant_K_module],
    :Air => [environmental_module, Ψ_air_module, constant_K_module],
)

# Rev her up

system = PlantModules.generate_system(struct_connections, func_connections, module_coupling, checkunits = false)
sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
@time sol = solve(prob);

# Plotting

plotgraph(sol, graphs[1], func_varname = :W)
plotgraph(sol, graphs[1], func_varname = :W, ylims = (0, 1))
plotgraph(sol, graphs[2], func_varname = :W)

plotgraph(sol, graphs[1], func_varname = :Ψ, struct_module = :Stem, ylims = (-0.15, -0.05))

plotgraph(sol, graphs[1], func_varname = :K_surface, struct_module = :Leaf)
plotgraph(sol, graphs[2], func_varname = :W)