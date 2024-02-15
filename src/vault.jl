# For functions I don't need anymore (but what if I end up needing them later ?!) #

# Apply a function to a node and all its descendants
function iteratedescendants(node, graph, func::Function; kwargs...)
	func(node; kwargs...)
	for chnode in PlantModules.children(node, graph)
		iteratedescendants(chnode, graph, func; kwargs...)
	end
end

# Default behaviour: start from graph root
iteratedescendants(graph, func::Function; kwargs...) = iteratedescendants(PlantModules.root(graph), graph, func; kwargs...)

# Fill in MTK_systems, MTK_connections, MTK_u0 for the given node
function add_MTK_info!(node; model_defaults, module_defaults, module_coupling, struct_connections,
	func_connections, MTK_systems, MTK_connections, MTK_u0)

	MTKsystem = getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params)
	# MTKconnections = getMTKconnections(node) # needs all connected nodes to be defined as MTKsystems already
end

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

root(graph::PlantGraphs.StaticGraph) = graph[graph.root]
root(graph::PlantGraphs.GraphNode) = graph

children(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = [graph[child_id] for child_id in node.children_id]