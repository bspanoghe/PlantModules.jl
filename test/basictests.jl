mutable struct Foo <: Node
    bar::Int
end

mutable struct Oof <: Node
    rab::String
end

function qux(; name, bar, baz, quux, corge, xyzzy)
    @parameters bar = bar baz = baz quux = quux corge = corge
    @variables t xyzzy(t) = xyzzy
    eqs = [
        bar + baz ~ quux * xyzzy + corge
    ]
    return ODESystem(eqs, t; name)
end

function snoo(; name, snee)
    @parameters snee = snee
    @variables t xyzzy1(t) xyzzy2(t) snaw(t)
    eqs = [
        snaw ~ xyzzy1 * xyzzy2 + snee
    ]
    return ODESystem(eqs, t; name)
end

PlantModules.attributes(node::Foo) = Dict([:bar => node.bar])
PlantModules.structmod(::Foo) = :Foo
PlantModules.structmod(::Oof) = :Oof
PlantModules.id(::Foo) = 1

func_module = qux
structmodule = :Foo
func_connections = [() => snoo]
module_coupling = [qux => [:Foo]]
default_params = (qux = (bar = 100, baz = 42, quux = 42, corge = 100), snoo = (snee = 100,))
default_u0s = (qux = (xyzzy = 100,),)
module_defaults = (Foo = (baz = 10, xyzzy = 10),)

## alter_defaults

default_changes = (baz = 100, quux = 100, xyzzy = 500)
altered_params, altered_u0s = alter_defaults(default_changes, default_params = default_params, default_u0s = default_u0s)
@test altered_params == (qux = (bar = 100, baz = 100, quux = 100, corge = 100), snoo = (snee = 100,))
@test altered_u0s == (qux = (xyzzy = 500,),)

## getnodeparamu0s

node = Foo(1)

@test PlantModules.getnodeparamu0s(node, structmodule, func_module, module_defaults, default_params) ==
    Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100])
@test PlantModules.getnodeparamu0s(node, structmodule, func_module, module_defaults, default_u0s) ==
    Dict([:xyzzy => 10])

## getMTKsystem

node = Foo(1)
testsystem = PlantModules.getMTKsystem(node, module_coupling, module_defaults, default_params, default_u0s)
@test get_name(testsystem) == :(Foo1)
@test sort(Symbol.(keys(get_defaults(testsystem)))) == sort(collect(keys(Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100, Symbol("xyzzy(t)") => 10]))))
@test sort(collect(values(get_defaults(testsystem)))) == sort(collect(values(Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100, Symbol("xyzzy(t)") => 10]))))

## get_func_connection
node = Foo(1)
nb_node = Foo(2)
testsystem = PlantModules.get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
@test get_name(testsystem) == :Foo1_Foo1
@test sort(Symbol.(keys(get_defaults(testsystem)))) == sort(collect(keys(Dict([:snee => 100]))))
@test sort(collect(values(get_defaults(testsystem)))) == sort(collect(values(Dict([:snee => 100]))))

## get_connection_info
graphnr = 1
graph = Foo(123) + (Foo(22), Foo(300))
node = PlantModules.nodes(graph)[1]
nb_nodes = PlantModules.nodes(graph)[2:3]
nb_node_graphnrs = [1, 1]
MTK_system_dicts = [Dict([PlantModules.nodes(graph)[nodenr].self_id =>
    qux(name = Symbol("Foo$nodenr"), bar = nodenr, baz = nodenr, quux = nodenr, corge = nodenr, xyzzy = nodenr)
    for nodenr in 1:3]
)]
get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
    [connection_MTK.xyzzy1 ~ node_MTK.xyzzy for connection_MTK in connection_MTKs]..., 
    [connection_MTK.xyzzy2 ~ nb_node_MTKs.xyzzy for (connection_MTK, nb_node_MTKs) in zip(connection_MTKs, nb_node_MTKs)]..., 
]

connection_MTKs, connection_equations = PlantModules.get_connection_info(node, graphnr, nb_nodes,
    nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts
)

@test connection_MTKs[1] == snoo(; name = Symbol("Foo" * string(PlantModules.id(node)) * "_" *
"Foo" * string(PlantModules.id(nb_nodes[1]))), snee = 100)
@test connection_MTKs[2] == snoo(; name = Symbol("Foo" * string(PlantModules.id(node)) * "_" *
"Foo" * string(PlantModules.id(nb_nodes[2]))), snee = 100)

@test connection_equations == get_connection_eqs(MTK_system_dicts[1][PlantModules.id(node)],
    [MTK_system_dicts[1][PlantModules.id(nb_nodes[1])], MTK_system_dicts[1][PlantModules.id(nb_nodes[2])]],
    connection_MTKs
)

## get_intergraph_neighbours
graphs = [Oof("1") + (Foo(2), Foo(3), Foo(4)), Oof("!")]
node_graphnr = 2
node = PlantModules.nodes(graphs[node_graphnr])[1]

### 1
intergraph_connections = [[1, 2] => (:Foo, :Oof)]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
@test ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
@test ig_neighbour_graphnrs == [1, 1, 1]

### 2
intergraph_connections = [[2, 1] => (:Foo, :Oof)]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
@test isempty(ig_neighbours)
@test isempty(ig_neighbour_graphnrs)

### 3
intergraph_connections = [[2, 1] => (:Oof, :Foo)]
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

intergraph_connections = [[1, 2] => (node1, node2) -> structmod(node1) == :Foo && attributes(node1)[:bar] in [2, 3] && structmod(node2) == :Oof]
ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, 2, graphs, intergraph_connections)
@test issetequal(ig_neighbours, [node2, node3])
@test ig_neighbour_graphnrs == [1, 1]

## get_nb_nodes
graphs = [Oof("1") + (Foo(2), Foo(3), Foo(4)), Foo(90)]
graphnr = 1
node = PlantModules.nodes(graphs[graphnr])[1]
intergraph_connections = [[1, 2] => (:Oof, :Foo)]

nb_nodes = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
@test nb_nodes[1] == vcat(PlantModules.nodes(graphs[1])[2:4], PlantModules.nodes(graphs[2])[1])
@test nb_nodes[2] == [1, 1, 1, 2]