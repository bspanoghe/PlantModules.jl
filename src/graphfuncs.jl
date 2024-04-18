"""
	nodes(graph)

Returns a `Vector` containing all nodes of the given graph.
"""
nodes(graph) = error("Function not yet defined for input type $(typeof(graph)).")

"""
	neighbours(node)

Returns a `Vector` containing the neighbours of the given node.
"""
neighbours(node) = error("Function not yet defined for input type $(typeof(node)).")

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

# Implementation for Base Julia representation as Dicts #

nodes(graph::Dict) = collect(values(graph))
neighbours(node::Dict, graph::Dict) = [graph[nbidx] for nbidx in node[:nb]]
attributes(node::Dict) = node[:at]
structmod(node::Dict) = node[:sm]
id(node::Dict) = node[:id]

# Implementation for PlantGraphs #

nodes(graph::PlantGraphs.StaticGraph) = graph.nodes.vals
nodes(node::PlantGraphs.Node) = [node] # Sometimes a graph is just a single node

neighbours(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = [graph[nb_id] for nb_id in (ismissing(node.parent_id) ? node.children_id : vcat(node.children_id.dict.keys, node.parent_id))]
neighbours(::PlantGraphs.Node, ::PlantGraphs.Node) = [] # node in a graph consisting of one node has no neighbours

function attributes(node::PlantGraphs.GraphNode)
	fields = fieldnames(typeof(node.data))
	
	if isempty(fields)
		return Dict([])
	end

	return Dict([field => getfield(node.data, field) for field in fields])
end

function attributes(node::PlantGraphs.Node)
	fields = fieldnames(typeof(node))
	
	if isempty(fields)
		return Dict([])
	end

	return Dict([field => getfield(node, field) for field in fields])
end

structmod(node::PlantGraphs.GraphNode) = typeof(node.data).name.name
structmod(node::PlantGraphs.Node) = typeof(node).name.name

id(node::PlantGraphs.GraphNode) = node.self_id
id(node::PlantGraphs.Node) = 1 # The entire graph only consists of one node if the input type is Node