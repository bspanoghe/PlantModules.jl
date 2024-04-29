# Required functions for all graphs #

"""
	nodes(graph)

Returns a `Vector` containing all nodes of the given graph.
"""
nodes(graph) = error("Function not yet defined for input type $(typeof(graph)).")

"""
	neighbours(node)

Returns a `Vector` containing the neighbours of the given node.
"""
neighbours(node, graph) = error("Function not yet defined for input types $(typeof(node)) and $(typeof(graph)).")

"""
	attributes(node)

Returns a `Dict` containing the functional variable names and values of the given node.
"""
attributes(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	structmod(node)

Returns the structural module corresponding with the node as a `Symbol`.
"""
structmod(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	id(node)

Returns the id of the node as an `Int`.
"""
id(node) = error("Function not yet defined for input type $(typeof(node)).")

## Implementation for Base Julia representation as Dicts

nodes(graph::Dict) = collect(values(graph))
neighbours(node::Dict, graph::Dict) = [graph[nbidx] for nbidx in (node[:parent_id] == -1 ? node[:child_ids] : vcat(node[:child_ids], node[:parent_id]))]
attributes(node::Dict) = node[:attributes]
structmod(node::Dict) = node[:type]
id(node::Dict) = node[:id]

## Implementation for PlantGraphs.jl

### GraphNode / StaticGraph (nodes combined into graph)

nodes(graph::PlantGraphs.StaticGraph) = graph.nodes.vals
neighbours(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = 
	[
		graph[nb_id]
	for nb_id in (ismissing(node.parent_id) ? node.children_id :
		vcat(node.children_id.dict.keys, node.parent_id))
]
function attributes(node::PlantGraphs.GraphNode)
	fields = fieldnames(typeof(node.data))
	
	if isempty(fields)
		return Dict([])
	end

	return Dict([field => getfield(node.data, field) for field in fields])
end
structmod(node::PlantGraphs.GraphNode) = typeof(node.data).name.name
id(node::PlantGraphs.GraphNode) = node.self_id

### Node (the graph is a single node)

nodes(node::PlantGraphs.Node) = [node]
neighbours(::PlantGraphs.Node, ::PlantGraphs.Node) = [] # node in a graph consisting of one node has no neighbours
function attributes(node::PlantGraphs.Node)
	fields = fieldnames(typeof(node))
	
	if isempty(fields)
		return Dict([])
	end

	return Dict([field => getfield(node, field) for field in fields])
end
endstructmod(node::PlantGraphs.Node) = typeof(node).name.name
id(node::PlantGraphs.Node) = 1 # The entire graph only consists of one node if the input type is Node

### MyPGNode (own implementation acquired by graph conversion)

nodes(node::MyPGNode) = [node]
neighbours(::MyPGNode, ::MyPGNode) = []
attributes(node::MyPGNode) = node[:attributes]
structmod(node::MyPGNode) = typeof(node).parameters[1]
id(node::MyPGNode) = 1


## Implementation for MultiScaleTreeGraph.jl

nodes(graph::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.descendants(graph, self = true)
neighbours(node::MultiScaleTreeGraph.Node, _) = isnothing(MultiScaleTreeGraph.parent(node)) ?
	MultiScaleTreeGraph.children(node) : vcat(MultiScaleTreeGraph.parent(node), MultiScaleTreeGraph.children(node))
attributes(node::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.node_attributes(node)
structmod(node::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.node_mtg(node).symbol
id(node::MultiScaleTreeGraph.Node) = MultiScaleTreeGraph.id(node)

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

root(graph::Dict) = [node for node in nodes(graph) if node[:parent_id] == -1][1]
children(node::Dict, graph::Dict) = [graph[child_id] for child_id in node[:child_ids]]
parent(node::Dict, graph::Dict) = node[:parent_id] == -1 ? nothing : graph[node[:parent_id]]

## Implementation for PlantGraphs.jl

## Implementation for MultiScaleTreeGraph.jl
