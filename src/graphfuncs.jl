"""
	nodes(graph)

Returns all nodes of the given graph.
"""
nodes(graph) = error("Function not yet defined for input type $(typeof(graph)).")

"""
	neighbours(node)

Returns the neighbours of the given node.
"""
neighbours(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	attributes(node)

Returns the attributes of the given node.
"""
attributes(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	nodetype(node)

Returns the type of the node, corresponding with the name of a structural module.
"""
nodetype(node) = error("Function not yet defined for input type $(typeof(node)).")

"""
	id(node)

Returns the id of the node.
"""
id(node) = error("Function not yet defined for input type $(typeof(node)).")

# Implementation for PlantGraphs #

nodes(graph::PlantGraphs.StaticGraph) = graph.nodes.vals
nodes(node::PlantGraphs.Node) = [node] # Sometimes a graph is just a single node

neighbours(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = [graph[nb_id] for nb_id in (ismissing(node.parent_id) ? node.children_id : vcat(node.children_id.dict.keys, node.parent_id))]

function attributes(node::PlantGraphs.GraphNode)
	fields = fieldnames(typeof(node.data))
	
	if isempty(fields)
		return Dict([])
	end

	fieldvalues = getfield.([node.data], fields...)
	return Dict([field => fieldvalue for (field, fieldvalue) in zip(fields, fieldvalues)])
end

function attributes(node::PlantGraphs.Node)
	fields = fieldnames(typeof(node))
	
	if isempty(fields)
		return Dict([])
	end

	fieldvalues = getfield.([node], fields...)
	return Dict([field => fieldvalue for (field, fieldvalue) in zip(fields, fieldvalues)])
end

nodetype(node::PlantGraphs.GraphNode) = string(node.data) |> x -> split(x, '(')[1] |> Symbol
nodetype(node::PlantGraphs.Node) = string(node) |> x -> split(x, '(')[1] |> Symbol

id(node::PlantGraphs.GraphNode) = node.self_id
id(node::PlantGraphs.Node) = 1 # The entire graph only consists of one node if the input type is Node