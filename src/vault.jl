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
0
## iteratedescendants 
testgraph = Foo(1) + (Foo(2), Foo(3) + (Foo(7), Foo(9)))
bars = Int[]
iteratedescendants(testgraph, (x; extra) -> push!(bars, x.data.bar + extra), extra = 3)
bars == [1, 2, 3, 7, 9] .+ 3


# get MTK systems and the vector of equations defining connections between MTK systems of all nodes between two graphs as specified in the intergraph connections
function get_intergraph_connection_info(graphs, intergraph_connection, func_connections,
	get_connection_eqs, default_params, default_u0s, MTK_system_dicts)

	graph1, graph2 = graphs[intergraph_connection[1]]
	MTK_system_dict, nb_MTK_system_dict = MTK_system_dicts[intergraph_connection[1]]

	nodes = [node for node in PlantModules.nodes(graph1) if PlantModules.nodetype(node) in intergraph_connection[2]]
	nb_nodes = [node for node in PlantModules.nodes(graph2) if PlantModules.nodetype(node) in intergraph_connection[2]]
	
	if isempty(nodes) || isempty(nb_nodes) # no neighbours no connection info
		return ODESystem[], Equation[]
	end

	full_connection_MTKs = Vector{ODESystem}(undef, length(nodes) * length(nb_nodes)) # we're fully connecting the valid nodes from both graphs
	full_connection_equations = Equation[]

	for (node_nr, node) in enumerate(nodes)
		node_MTK = MTK_system_dict[PlantModules.id(node)]
		
		nb_node_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
		connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes))

		for (nb_nr, nb_node) in enumerate(nb_nodes)
			nb_node_MTKs[nb_nr] = nb_MTK_system_dict[PlantModules.id(nb_node)]
			connection_MTKs[nb_nr] = get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
		end
	
		push!(full_connection_equations, get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs)...)
		full_connection_MTKs[node_nr:node_nr+length(nb_nodes)-1] = connection_MTKs
	end

	return full_connection_MTKs, full_connection_equations
end

intergraph_connections = struct_connections[2]
for intergraph_connection in intergraph_connections
	intergraph_connection_MTKs, intergraph_connection_equations = get_intergraph_connection_info(graphs, intergraph_connection, func_connections,
		get_connection_eqs, default_params, default_u0s, MTK_system_dicts
	)
	push!(connection_MTKs, intergraph_connection_MTKs...)
	push!(connection_equations, intergraph_connection_equations...)
end