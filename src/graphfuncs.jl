"""
	root(graph)

Returns the root node of the given graph.
"""
root(graph) = error("Function not yet defined for given input type.")
root(graph::PlantGraphs.StaticGraph) = graph[1]
root(graph::PlantGraphs.GraphNode) = graph

"""
	children(node)

Returns the children nodes of the given node.
"""
children(node) = error("Function not yet defined for given input type.")
children(node::PlantGraphs.GraphNode) = PlantGraphs.children(node)
