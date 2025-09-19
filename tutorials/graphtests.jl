using Pkg; Pkg.activate("./tutorials")
using Graphs, Plots, GraphRecipes

begin
    a = SimpleGraph(10)
    for i in 1:10
        add_edge!(a, i-1, i)
    end
end
a

neighbors(a, 1)
graphplot(a)

[i for i in BFSIterator(a, 1)]

vertices(a)
edges(a) |> collect
edgetype(a)
has_edge(a, 1, 2)
has_vertex(a, 1)
@which inneighbors(a, 1)
@btime outneighbors(a, 1)
ne(a)
nv(a)
is_directed(a)


import Graphs: edges, edgetype, has_edge, has_vertex, inneighbors, ne, nv, outneighbors, vertices, is_directed

struct PMVertex
    id::Int64
    type::Symbol
    attributes::Dict
end

id(v::PMVertex) = v.id
type(v::PMVertex) = v.type
attributes(v::PMVertex) = v.attributes

struct PMGraph{T} <: AbstractGraph{T}
    vertices::Vector{T}
    neighbors::Dict{T, Vector{T}}
    ne::Integer
    vertexdict::Dict{T, PMVertex}
end

vertices(g::PMGraph) = g.vertices
edges(::PMGraph) = nothing
ne(g::PMGraph) = g.ne
vertexdict(g::PMGraph{T}, v::T) where {T} = g.vertexdict[v]

edgetype(g::PMGraph) = eltype(edges(g))
neighbors(g::PMGraph, v::PMVertex) = g.neighbors[v]
inneighbors(g::PMGraph, v::PMVertex) = neighbors(g, v)
outneighbors(g::PMGraph, v::PMVertex) = neighbors(g, v)
has_edge(g::PMGraph, s::PMVertex, d::PMVertex) = d in neighbors(g, s)
has_vertex(g::PMGraph, v::PMVertex) = v in vertices(g)
nv(g::PMGraph) = length(vertices(g))
is_directed(g::PMGraph) = false
is_directed(PMGraph) = false

getnodes(g::PMGraph) = [vertexdict(g, v) for v in vertices(g)]
getneighbours(v::PMVertex, g::PMGraph) = [vertexdict(g, v) for v in neighbors(g, id(v))]
getattributes(v::PMVertex) = attributes(v)
getid(v::PMVertex) = id(v)
getstructmod(v::PMVertex) = type(v)



vertexlist = [PMVertex(i, :horse, Dict(:cheese => rand())) for i in 1:10]
function addneighbours!(nd::Dict{T, Vector{T}}, v1::T, v2::T) where {T}
    v1_nbs = get(nd, v1, T[])
    v2 in values(v1_nbs) && return false
    v2_nbs = get(nd, v2, T[])
    nd[v1] = [v1_nbs; v2]
    nd[v2] = [v2_nbs; v1]
    return true
end

begin
    neighbordict = Dict{PMVertex, Vector{PMVertex}}()
    num_edges = 0
    for _ in 1:10
        addedsomething = addneighbours!(neighbordict, vertexlist[ceil(Int64, 10*rand())], vertexlist[ceil(Int64, 10*rand())])
        addedsomething && (num_edges += 1)
    end
end

gg = PMGraph(vertexlist, neighbordict, num_edges)

vertices(gg)
degree(gg)

graphplot(gg)