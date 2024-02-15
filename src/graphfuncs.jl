"""
	nodes(graph)

Returns all nodes of the given graph.
"""
nodes(graph) = error("Function not yet defined for given input type.")

"""
	neighbours(node)

Returns the neighbours of the given node.
"""
neighbours(node) = error("Function not yet defined for given input type.")

"""
	attributes(node)

Returns the attributes of the given node.
"""
attributes(node) = error("Function not yet defined for given input type.")

"""
	nodetype(node)

Returns the type of the node, corresponding with the name of a structural module.
"""
nodetype(node) = error("Function not yet defined for given input type.")

"""
	id(node)

Returns the id of the node.
"""
id(node) = error("Function not yet defined for given input type.")

# Implementation for PlantGraphs #

function attributes(node::PlantGraphs.GraphNode)
	fields = fieldnames(typeof(node.data))
	
	if isempty(fields)
		return Dict([])
	end

	fieldvalues = getfield.([node.data], fields...)
	return Dict([field => fieldvalue for (field, fieldvalue) in zip(fields, fieldvalues)])
end

nodetype(node::PlantGraphs.GraphNode) = string(node.data) |> x -> split(x, '(')[1] |> Symbol

id(node::PlantGraphs.GraphNode) = node.self_id