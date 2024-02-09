"""
	root(graph)

Returns the root node of the given graph.
"""
root(graph) = error("Function not yet defined for given input type.")


"""
	children(node)

Returns the children nodes of the given node.
"""
children(node) = error("Function not yet defined for given input type.")

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

# Implementation for PlantGraphs #

root(graph::PlantGraphs.StaticGraph) = graph[graph.root]
root(graph::PlantGraphs.GraphNode) = graph

children(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = [graph[child_id] for child_id in node.children_id]

function attributes(node::PlantGraphs.GraphNode)
	fields = fieldnames(typeof(node.data))
	fieldvalues = getfield.([node.data], fields...)
	return [field => fieldvalue for (field, fieldvalue) in zip(fields, fieldvalues)]
end

nodetype(node::PlantGraphs.GraphNode) = string(node.data) |> x -> split(x, '(')[1] |> Symbol