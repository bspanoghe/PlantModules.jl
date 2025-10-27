# Structure definition

Base.@kwdef mutable struct Forest <: PlantGraphs.Node
    N::Int = 20
    P::Int = 10
    δ::Float64 = 1.0
end

Base.@kwdef mutable struct Grassland <: PlantGraphs.Node
    N::Int = 100
    P::Int = 5
end

Base.@kwdef mutable struct Cave <: PlantGraphs.Node end

graph = Grassland() + Grassland(N = 50, P = 3) + (Forest(δ = 1.8), Forest(N = 5, P = 0, δ = 2.5))
graph2 = Cave() + Forest() + Grassland()
graphs = [graph, graph2]
intergraph_connections = [[1, 2] => (:Forest, :Cave)]
plantstructure = PlantStructure(graphs, intergraph_connections)

# PlantModules graph functions

node1, node2, node3, node4 = getnodes(graph)

## getnodes
@test issetequal(PlantModules.getnodes(graph), [node1, node2, node3, node4])

## getneighbors
@test issetequal(PlantModules.getneighbors(node2, graph), [node1, node3, node4])

## getattributes
@test issetequal(PlantModules.getattributes(node4), [:P => 0, :N => 5, :δ => 2.5])

## getstructmod
@test PlantModules.getstructmod(node1) == :Grassland

## getid
@test allunique(PlantModules.getid.([node1, node2, node3, node4]))

# PlantStructure functions

## get_nb_nodes
node = node4
graphnr = 1
nb_nodes = PlantModules.get_nb_nodes(node, graphnr, graphs, intergraph_connections)
@test issetequal(nb_nodes, [node2, PlantModules.getnodes(graphs[2])[1]])