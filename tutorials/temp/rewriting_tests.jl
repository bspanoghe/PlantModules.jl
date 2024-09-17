using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs
# using MultiScaleTreeGraph
using ModelingToolkit, DifferentialEquations, Unitful
using Plots; import GLMakie.draw
using BenchmarkTools

# Plant structure

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
    leaf_D = [2.0, 1.0, 0.1]
    parent_D = data(x).D
    child_D = [parent_D[1] / sqrt(2), 1.0]
    return Stem(parent_D) + (Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)), Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)))
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

plant = Graph(axiom = axiom, rules = (rule1,))
num_iterations = 2
for _ in 1:num_iterations
    rewrite!(plant)
	growify!(plant, 1.1)
end

draw(plant)
convert_to_MTG(plant) |> PlantModules.MultiScaleTreeGraph.DataFrame

C_root = 300e-6
C_stem = 400e-6
C_leaf = 450e-6

module_defaults = (
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 15.0], ϕ_D = [0.8, 0.03]), D = [1.5, 10], M = C_stem),
	Leaf = (shape = PlantModules.Cuboid(ϵ_D = [10.0, 10.0, 0.1], ϕ_D = [0.5, 0.5, 0.01]), M = C_leaf),
	Soil = (W_max = 10000.0, T = 288.15),
	Air = (W_r = 0.8,)
)

module_coupling = [
	PlantModules.hydraulic_module => [:Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

soil_graph = Soil()
air_graph = Air()

graphs = [plant, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (PlantModules.root(plant), :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)]

struct_connections = [graphs, intergraph_connections]


function sizedep_hydraulic_connection(; name, K, connection_shape)
    @parameters (
        K_s = K_s, [description = "Specific hydraulic conductivity of connection", unit = u"g / hr / MPa / cm^2"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
		D_1(t), [description = "Dimensions of compartment 1", unit = u"cm^2"],
		D_2(t), [description = "Dimensions of compartment 2", unit = u"cm^2"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1)
		K ~ K_s * (PlantModules.surface_area(connection_shape, D_1) + PlantModules.surface_area(connection_shape, D_2))/2
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order) = [ 
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
		connection_MTK.D_1 ~ node_MTK.D,
        connection_MTK.D_2 ~ nb_node_MTK.D,
    ]

    return ODESystem(eqs, t; name), get_connection_eqset
end

connecting_modules = [
	(:Soil, :Stem) => (PlantModules.hydraulic_connection, [:K => 80]),
	(:Stem, :Stem) => (PlantModules.hydraulic_connection, [:K => 80]),
	(:Stem, :Leaf) => (PlantModules.hydraulic_connection, [:K => 40]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-3]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-3]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

func_connections = [connecting_modules, PlantModules.multi_connection_eqs]

system = PlantModules.generate_system(PlantModules.default_params, PlantModules.default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
@time sol = solve(prob);

# PlantModules.plotgraph(sol, graphs[1], func_varname = :W)
# ylims!(0, 1)
# PlantModules.plotgraph(sol, graphs[2], func_varname = :W)
# PlantModules.plotgraph(sol, graphs[1], func_varname = :P)
# PlantModules.plotgraph(sol, graphs[1], func_varname = :Ψ)
# PlantModules.plotgraph(sol, graphs[1], func_varname = :ΣF)

# [sol[Symbol(string(sys.name) * "₊W")] for sys in sol.prob.f.sys.systems if !occursin("_", string(sys.name))] .|> minimum |> minimum
# [sol[Symbol(string(sys.name) * "₊ΣF")] for sys in sol.prob.f.sys.systems if !occursin("_", string(sys.name))] .|> maximum |> maximum