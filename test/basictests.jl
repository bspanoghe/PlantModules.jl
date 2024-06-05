@variables t
d = Differential(t)

mutable struct Stem <: Node
    age::Int
end

mutable struct Leaf <: Node
    color::String
end

function water_module(; name, age, psi_m, psi_p, psi_s, water)
    @parameters age = age psi_m = psi_m psi_p = psi_p psi_s = psi_s
    @variables water(t) = water
    eqs = [
        d(water) ~ (psi_m + psi_p + psi_s) * age
    ]
    return ODESystem(eqs, t; name)
end

function carbon_module(; name, ps_rate)
    @parameters ps_rate = ps_rate
    @variables co2(t) metabolite(t)
    eqs = [
        d(metabolite) ~ ps_rate*co2
    ]
    return ODESystem(eqs, t; name)
end

function co2_module(; name, co2)
    @variables co2(t) = co2
    eqs = [
        d(co2) ~ sin(t)
    ]
    return ODESystem(eqs, t; name)
end

# PlantModules.attributes(node::Stem) = Dict([:age => node.age])
# PlantModules.structmod(::Stem) = :Stem
# PlantModules.structmod(::Leaf) = :Leaf
# PlantModules.id(::Stem) = 1

func_module = water_module
structmodule = :Stem
func_connections = [() => carbon_module]
module_coupling = [water_module => [:Stem]]
default_params = (water_module = (age = 100, psi_m = 42, psi_p = 42, psi_s = 100), carbon_module = (ps_rate = 100,))
default_u0s = (water_module = (water = 100,),)
module_defaults = (Stem = (psi_m = 10, water = 10),)

## alter_defaults

default_changes = (psi_m = 100, psi_p = 100, water = 500)
altered_params, altered_u0s = alter_defaults(default_changes, default_params = default_params, default_u0s = default_u0s)
@test altered_params == (water_module = (age = 100, psi_m = 100, psi_p = 100, psi_s = 100), carbon_module = (ps_rate = 100,))
@test altered_u0s == (water_module = (water = 500,),)

## getnodeparamu0s

node = Stem(1)

@test PlantModules.getnodeparamu0s(node, structmodule, func_module, module_defaults, default_params) ==
    Dict([:age => 1, :psi_m => 10, :psi_p => 42, :psi_s => 100])
@test PlantModules.getnodeparamu0s(node, structmodule, func_module, module_defaults, default_u0s) ==
    Dict([:water => 10])

## collapse
checkunits = false

sys1 = co2_module(name = :sys1, co2 = 10)
sys2 = carbon_module(name = :sys2, ps_rate = 5)
sys12 = PlantModules.collapse([sys1, sys2], name = :sys12, checkunits = checkunits)
@test Symbol(get_iv(sys12)) == Symbol(get_iv(sys1))
@test issetequal(get_eqs(sys12), [get_eqs(sys1)..., get_eqs(sys2)...])
@test issetequal(unknowns(sys12), [unknowns(sys1)..., unknowns(sys2)...])
@test issetequal(parameters(sys12), [parameters(sys1)..., parameters(sys2)...])
@test issetequal(get_defaults(sys12), [get_defaults(sys1)..., get_defaults(sys2)...])

## getMTKsystem

node = Stem(1)
checkunits = false
testsystem = PlantModules.getMTKsystem(node, module_coupling, module_defaults, default_params, default_u0s, checkunits)
@test get_name(testsystem) == :(Stem1)
@test sort(Symbol.(keys(get_defaults(testsystem)))) ==
    sort(collect(keys(Dict([:age => 1, :psi_m => 10, :psi_p => 42, :psi_s => 100, Symbol("water(t)") => 10]))))
@test sort(collect(values(get_defaults(testsystem)))) ==
    sort(collect(values(Dict([:age => 1, :psi_m => 10, :psi_p => 42, :psi_s => 100, Symbol("water(t)") => 10]))))

## get_func_connection
node = Stem(1)
nb_node = Stem(2)
testsystem = PlantModules.get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
@test get_name(testsystem) == :Stem1_Stem1
@test sort(Symbol.(keys(get_defaults(testsystem)))) == sort(collect(keys(Dict([:ps_rate => 100]))))
@test sort(collect(values(get_defaults(testsystem)))) == sort(collect(values(Dict([:ps_rate => 100]))))

## get_connection_info
graphnr = 1
graph = Stem(123) + (Stem(22), Stem(300))
node = PlantModules.nodes(graph)[1]
nb_nodes = PlantModules.nodes(graph)[2:3]
nb_node_graphnrs = [1, 1]
MTK_system_dicts = [Dict([PlantModules.nodes(graph)[nodenr].self_id =>
    water_module(name = Symbol("Stem$nodenr"), age = nodenr, psi_m = nodenr, psi_p = nodenr, psi_s = nodenr, water = nodenr)
    for nodenr in 1:3]
)]
get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
    [connection_MTK.co2 ~ node_MTK.water for connection_MTK in connection_MTKs]..., 
    [connection_MTK.metabolite ~ nb_node_MTKs.water for (connection_MTK, nb_node_MTKs) in zip(connection_MTKs, nb_node_MTKs)]..., 
]

connection_MTKs, connection_equations = PlantModules.get_connection_info(node, graphnr, nb_nodes,
    nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts
)

@test connection_MTKs[1] == carbon_module(; name = Symbol("Stem" * string(PlantModules.id(node)) * "_" *
"Stem" * string(PlantModules.id(nb_nodes[1]))), ps_rate = 100)
@test connection_MTKs[2] == carbon_module(; name = Symbol("Stem" * string(PlantModules.id(node)) * "_" *
"Stem" * string(PlantModules.id(nb_nodes[2]))), ps_rate = 100)

@test connection_equations == get_connection_eqs(MTK_system_dicts[1][PlantModules.id(node)],
    [MTK_system_dicts[1][PlantModules.id(nb_nodes[1])], MTK_system_dicts[1][PlantModules.id(nb_nodes[2])]],
    connection_MTKs
)

## get_intergraph_neighbours
graphs = [Leaf("1") + (Stem(2), Stem(3), Stem(4)), Leaf("!")]
node_graphnr = 2
node = PlantModules.nodes(graphs[node_graphnr])[1]

### 1
intergraph_connections = [[1, 2] => (:Stem, :Leaf)]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
@test ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
@test ig_neighbour_graphnrs == [1, 1, 1]

### 2
intergraph_connections = [[2, 1] => (:Stem, :Leaf)]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
@test isempty(ig_neighbours)
@test isempty(ig_neighbour_graphnrs)

### 3
intergraph_connections = [[2, 1] => (:Leaf, :Stem)]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
@test ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
@test ig_neighbour_graphnrs == [1, 1, 1]

### 4

node2 = PlantModules.nodes(graphs[1])[2]
node3 = PlantModules.nodes(graphs[1])[3]
intergraph_connections = [[1, 2] => ([node2, node3], [node])]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, 2, graphs, intergraph_connections)
@test issetequal(ig_neighbours, [node2, node3])
@test ig_neighbour_graphnrs == [1, 1]

### 5

intergraph_connections = [[1, 2] => (node1, node2) -> structmod(node1) == :Stem && attributes(node1)[:age] in [2, 3] && structmod(node2) == :Leaf]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, 2, graphs, intergraph_connections)
@test issetequal(ig_neighbours, [node2, node3])
@test ig_neighbour_graphnrs == [1, 1]

## get_nb_nodes
graphs = [Leaf("1") + (Stem(2), Stem(3), Stem(4)), Stem(90)]
graphnr = 1
node = PlantModules.nodes(graphs[graphnr])[1]
intergraph_connections = [[1, 2] => (:Leaf, :Stem)]

nb_nodes = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
@test nb_nodes[1] == vcat(PlantModules.nodes(graphs[1])[2:4], PlantModules.nodes(graphs[2])[1])
@test nb_nodes[2] == [1, 1, 1, 2]