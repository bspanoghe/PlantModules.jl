# # Graphs.jl `AbstractGraph` implementation with numerical attributes

# ## Types & constructors

struct PMVertex{T}
    id::T
    type::Symbol
    attributes::Dict
end

id(v::PMVertex) = v.id
type(v::PMVertex) = v.type
attributes(v::PMVertex) = v.attributes

struct PMEdge{T} <: AbstractEdge{T}
    src::T
    dst::T
end

src(e::PMEdge) = e.src
dst(e::PMEdge) = e.dst

"""
	PlantStructure{T} <: AbstractGraph{T}

The internal graph implementation used for [`generate_system`](@ref).

# Constructor

	PlantStructure(graphs::Vector, intergraph_connections::Vector; return_id_conversions::Bool = false)

Connect graphs together into a single `PlantStructure` graph.

Warning: This function defines its own node ids to simplify the integration with the rest of the `Graphs.jl` ecosystem.

Inputs:
- `graphs::Vector`: A vector of graphs representing the plants and their environment.
- `intergraph_connections::Vector`: A vector of pairs between 2 indexes of above graphs and how they are linked.
- `return_id_conversions::Bool`: Should the function return a dictionary mapping the newly defined node ids to the original nodes?

# Alternatives

Other graph implementations can also be used by extending the following functions: 
- [`getnodes`](@ref)
- [`getneighbors`](@ref)
- [`getattributes`](@ref)
- [`getstructmod`](@ref)
- [`getid`](@ref)
"""
struct PlantStructure{T} <: AbstractGraph{T}
    vertices::Vector{T}
    edges::Vector{PMEdge{T}}
    pmvertexdict::Dict{T, PMVertex{T}}
    neighbordict::Dict{T, Vector{T}}
end

function PlantStructure(vertices::Vector{T}, edges::Vector{PMEdge{T}}, pmvertexdict::Dict{T, PMVertex{T}}) where {T}
    neighbordict = Dict{T, Vector{T}}()
    for e in edges
        v1, v2 = src(e), dst(e)
        v1_nbs = get(neighbordict, v1, T[])
        v2_nbs = get(neighbordict, v2, T[])
        v2 in v1_nbs || (neighbordict[v1] = [v1_nbs; v2])
        v1 in v2_nbs || (neighbordict[v2] = [v2_nbs; v1])
    end
    return PlantStructure(vertices, edges, pmvertexdict, neighbordict)
end

vertices(g::PlantStructure) = g.vertices
edges(g::PlantStructure) = g.edges
pmvertex(g::PlantStructure{T}, v::T) where {T} = g.pmvertexdict[v]
neighbors(g::PlantStructure{T}, v::T) where {T} = g.neighbordict[v]
neighbors(g::PlantStructure{T}, v::T) where {T <: Integer} = g.neighbordict[v]

edgetype(g::PlantStructure) = eltype(edges(g))
inneighbors(g::PlantStructure{T}, v::T) where {T} = neighbors(g, v)
outneighbors(g::PlantStructure{T}, v::T) where {T} = neighbors(g, v)
has_edge(g::PlantStructure{T}, s::T, d::T) where {T} = d in neighbors(g, s)
has_vertex(g::PlantStructure{T}, v::T) where {T} = v in vertices(g)
nv(g::PlantStructure) = length(vertices(g))
ne(g::PlantStructure) = length(edges(g))
is_directed(::PlantStructure) = false
is_directed(::Type{<:PlantStructure}) = false

# ## Construct from other graphs

function PlantStructure(graphs::Vector, intergraph_connections::Vector; return_id_conversions::Bool = false)
    T = graphs[1] |> getnodes |> x -> getid(x[1]) |> typeof
    vertices = T[]
    edges = PMEdge{T}[]
    pmvertexdict = Dict{T, PMVertex{T}}()

    allnodes = [node for graph in graphs for node in getnodes(graph)]
    id_dict = Pair.(allnodes, eachindex(allnodes)) |> Dict

    # connect nodes to single graph
    for (graphnr, graph) in enumerate(graphs)
        for node in getnodes(graph)
            node_id = id_dict[node]
            push!(vertices, node_id)
            pmvertexdict[node_id] = PMVertex(node_id, getstructmod(node), getattributes(node))

            nb_nodes = get_nb_nodes(node, graphnr, graphs, intergraph_connections)
            append!(edges, [PMEdge(node_id, id_dict[nb_node]) for nb_node in nb_nodes])
        end
    end

    return_id_conversions && return PlantStructure(vertices, edges, pmvertexdict), id_dict
    return PlantStructure(vertices, edges, pmvertexdict)
end

PlantStructure(graph; return_id_conversions::Bool = false) = PlantStructure([graph], []; return_id_conversions)

# get neighbouring nodes of a node both from the same graph and all connected graphs
function get_nb_nodes(node, graphnr, graphs, intergraph_connections)
    graph = graphs[graphnr]

    intra_nb_nodes = getneighbors(node, graph)
    inter_nb_nodes = get_intergraph_neighbours(node, graphnr, graphs, intergraph_connections)

    nb_nodes = vcat(intra_nb_nodes, inter_nb_nodes)
    # isempty(nb_nodes) && error("No neighbours found for node $node.") #!

    return nb_nodes
end

# get neighbouring nodes from different graphs
## go over all intergraph connections
function get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
    nb_nodes = []

    # Go over all graphs and add nodes to neighbouring nodes if the graphs are connected
    for intergraph_connection in intergraph_connections
        if node_graphnr in first(intergraph_connection)
            node_idx, nb_idx = first(intergraph_connection)[1] == node_graphnr ? (1, 2) : (2, 1)
            nb_graphnr = first(intergraph_connection)[nb_idx]
            nb_graph = graphs[nb_graphnr]
            connection = intergraph_connection[2]
            if connection isa Tuple
                _nb_nodes = _get_intergraph_neighbours(node, nb_graph, connection[node_idx], connection[nb_idx])
            else
                nb_first = nb_idx == 1
                _nb_nodes = _get_intergraph_neighbours(node, nb_graph, connection, nb_first)
            end

            append!(nb_nodes, _nb_nodes)
        end
    end

    return nb_nodes
end

## For a connection where the structural module or specific nodes are specified
function _get_intergraph_neighbours(node, nb_graph, node_connection, nb_connection)
    if connection_check(node, node_connection)
        nb_nodes = [nb_node for nb_node in getnodes(nb_graph) if connection_check(nb_node, nb_connection)]
    else
        nb_nodes = []
    end

    return nb_nodes
end

connection_check(node, connection) = (node == connection) # connection is a node
connection_check(node, connection::Symbol) = (getstructmod(node) == connection) # connection is a structural module
connection_check(node, connection::Vector) = node in connection # connection is a collection of nodes

## For a connection with a user-defined filter function
function _get_intergraph_neighbours(node, nb_graph, connection_func::Function, nb_first::Bool)
    # Main.@infiltrate getstructmod(node) == :Air
    if nb_first
        nb_nodes = [nb_node for nb_node in getnodes(nb_graph) if connection_func(nb_node, node)]
    else
        nb_nodes = [nb_node for nb_node in getnodes(nb_graph) if connection_func(node, nb_node)]
    end

    return nb_nodes
end

# # PlantGraphs.jl `Node` implementation for graph conversion

struct MyPGNode{T} <: PlantGraphs.Node
    attributes::Dict
end

MyPGNode(type::Symbol, attributes) = MyPGNode{type}(attributes)
