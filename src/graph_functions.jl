# Required functions for all graphs #

"""
	getnodes(graph)

Returns a `Vector` containing all nodes of the given graph.
"""
getnodes(graph) = error("Function not yet defined for input type $(typeof(graph)).")

"""
	getneighbours(node)

Returns a `Vector` containing the neighbours of the given node.
"""
getneighbours(node, graph) = error("Function not yet defined for input types $(typeof(node)) and $(typeof(graph)).")

"""
	getattributes(node)

Returns a `Dict` containing the functional variable names and values of the given node.
"""
getattributes(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	getstructmod(node)

Returns the structural module corresponding with the node as a `Symbol`.
"""
getstructmod(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	getid(node)

Returns the id of the node as an `Int`.
"""
getid(node) = error("Function not yet defined for input type $(typeof(node)).")

## Implementation for Base Julia representation as Dicts

getnodes(graph::Dict) = collect(values(graph))
getneighbours(node::Dict, graph::Dict) = [graph[nbidx] for nbidx in (node[:parent_id] == -1 ? node[:child_ids] : vcat(node[:child_ids], node[:parent_id]))]
getattributes(node::Dict) = node[:attributes]
getstructmod(node::Dict) = node[:type]
getid(node::Dict) = node[:id]

## Implementation for PlantGraphs.jl

### GraphNode / StaticGraph (nodes combined into graph)

getnodes(graph::PlantGraphs.StaticGraph) = values(graph.nodes) |> collect
getneighbours(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = 
	[
		graph[nb_id]
	for nb_id in (ismissing(node.parent_id) ? node.children_id :
		vcat(node.children_id.dict.keys, node.parent_id))
]
function getattributes(node::PlantGraphs.GraphNode)
	fields = fieldnames(typeof(node.data))
	
	if isempty(fields)
		return Dict([])
	end

	return Dict([field => getfield(node.data, field) for field in fields])
end
getstructmod(node::PlantGraphs.GraphNode) = typeof(node.data).name.name
getid(node::PlantGraphs.GraphNode) = node.self_id


### Graph (dynamic graphs used in rewriting)

getnodes(graph::PlantGraphs.Graph) = getnodes(graph.graph)
getneighbours(node::PlantGraphs.GraphNode, graph::PlantGraphs.Graph) = getneighbours(node, graph.graph)

### Node (the graph is a single node)

getnodes(node::PlantGraphs.Node) = [node]
getneighbours(::PlantGraphs.Node, ::PlantGraphs.Node) = [] # node in a graph consisting of one node has no neighbours
function getattributes(node::PlantGraphs.Node)
	fields = fieldnames(typeof(node))
	
	if isempty(fields)
		return Dict([])
	end

	return Dict([field => getfield(node, field) for field in fields])
end
getstructmod(node::PlantGraphs.Node) = typeof(node).name.name
getid(node::PlantGraphs.Node) = 1 # The entire graph only consists of one node if the input type is Node

## MyPGNode (own implementation acquired by graph conversion)

getnodes(node::MyPGNode) = [node]
getneighbours(::MyPGNode, ::MyPGNode) = []
getattributes(node::MyPGNode) = node[:attributes]
getstructmod(node::MyPGNode) = typeof(node).parameters[1]
getid(node::MyPGNode) = 1


## Implementation for MultiScaleTreeGraph.jl

getnodes(graph::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.descendants(graph, self = true)
getneighbours(node::MultiScaleTreeGraph.Node, _) = isnothing(MultiScaleTreeGraph.parent(node)) ?
	MultiScaleTreeGraph.children(node) : vcat(MultiScaleTreeGraph.parent(node), MultiScaleTreeGraph.children(node))
getattributes(node::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.node_getattributes(node)
getstructmod(node::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.node_mtg(node).symbol |> Symbol
getid(node::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.node_id(node)

# Functions for Tree Graphs #

"""
	root(graph)

Returns the root node of the graph.
"""
root(graph) = error("Function not yet defined for input type $(typeof(graph)).")

"""
	children(node, graph)

Returns the child nodes of a node.
"""
children(node, graph) = error("Function not yet defined for input types $(typeof(node)) and $(typeof(graph)).")

"""
	parent(node, graph)

Returns the parent node of a node.
"""
parent(node, graph) = error("Function not yet defined for input types $(typeof(node)) and $(typeof(graph)).")

## Implementation for Base Julia representation as Dicts

root(graph::Dict) = [node for node in getnodes(graph) if node[:parent_id] == -1][1]
children(node::Dict, graph::Dict) = [graph[child_id] for child_id in node[:child_ids]]
parent(node::Dict, graph::Dict) = node[:parent_id] == -1 ? nothing : graph[node[:parent_id]]

## Implementation for PlantGraphs.jl

### StaticGraph

root(graph::PlantGraphs.StaticGraph) = graph[graph.root]
children(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = [graph[child_id] for child_id in node.children_id]
parent(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = graph[node.parent_id]

### (dynamic) Graph

root(graph::PlantGraphs.Graph) = root(graph.graph)
children(node::PlantGraphs.GraphNode, graph::PlantGraphs.Graph) = children(node, graph.graph)
parent(node::PlantGraphs.GraphNode, graph::PlantGraphs.Graph) = parent(node, graph.graph)

## Implementation for MultiScaleTreeGraph.jl

root(graph::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.get_root(graph)
children(node::MultiScaleTreeGraph.Node, _) = getfield(node, :children)
parent(node::MultiScaleTreeGraph.Node, _) = MultiScaleTreeGraph.parent(node)
