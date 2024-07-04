"""
	generate_system(default_params::NamedTuple, default_u0s::NamedTuple, module_defaults::NamedTuple,
		module_coupling::Vector, struct_connections::Vector, func_connections::Vector; checkunits::Bool)

Creates a MTK system based on a set of structural and functional modules and how those are connected.

# Arguments
- `default_params`: Model-wide default parameter values.
- `default_u0s`: Model-wide default initial values.
- `module_defaults`: Module-specific default values of both parameters and initial values.
- `module_coupling`: Coupling between functional and structural modules.
- `struct_connections`: One or more graphs specifying how the structural modules are connected, alongside the connections between graphs.
- `func_connections`: Additional functional information about the connections between structural modules.
- `checkunits`: Should the model check the units of the given equations? Defaults to `true`.
"""
function generate_system(default_params::NamedTuple, default_u0s::NamedTuple, module_defaults::NamedTuple,
	module_coupling::Vector, struct_connections::Vector, func_connections::Vector; checkunits::Bool = true
	)

	graphs, intergraph_connections = struct_connections
	connecting_modules, multi_connection_eqset = func_connections

	MTK_system_dicts = get_MTK_system_dicts(
		graphs, module_coupling, module_defaults, default_params, default_u0s, checkunits
	) # Vector of a dict per graph with (node id => node MTK system)
	MTK_systems = vcat(collect.(values.(MTK_system_dicts))...) # node MTK systems

	connection_MTKs = ODESystem[] # MTK systems of edges
	connection_eqsets = Equation[] # equations linking nodes with edges

	for (graphnr, graph) in enumerate(graphs) # go over all graphs
		for node in PlantModules.nodes(graph) # go over every node in the graph
			
			node_id = PlantModules.id(node)
			nb_nodes, nb_node_graphnrs = get_nb_nodes(node, graphs, graphnr, intergraph_connections) # collect neighbour nodes
			current_connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes)) # MTKs of current node's connections

			for (nb_idx, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs)) # go over all neighbours of the node
				connecting_module, reverse_order = get_connecting_module(node, nb_node, connecting_modules)
				connection_MTK, connection_eqset = get_connection_info( 
					node, graphnr, nb_node, nb_node_graphnr, connecting_module,
					reverse_order, default_params, default_u0s, MTK_system_dicts
				)
				current_connection_MTKs[nb_idx] = connection_MTK
				append!(connection_eqsets, connection_eqset)
			end
			append!(connection_MTKs, current_connection_MTKs)
			append!(connection_eqsets, multi_connection_eqset(MTK_system_dicts[graphnr][node_id], current_connection_MTKs))

		end
	end

	system = ODESystem(connection_eqsets, get_iv(MTK_systems[1]), name = :system,
		systems = vcat(MTK_systems, connection_MTKs), checks = checkunits
	) # combine all subsystems together with equations linking nodes with edges

	return system
end

# get node idx => node MTK system
function get_MTK_system_dicts(graphs, module_coupling, module_defaults, default_params, default_u0s, checkunits) 
	return [
		[PlantModules.id(node) => 
			getMTKsystem(node, module_coupling, module_defaults, default_params, default_u0s, checkunits)
			for node in PlantModules.nodes(graph)
		] |> Dict for graph in graphs
	]
end

# Get MTK system corresponding with node
function getMTKsystem(node, module_coupling, module_defaults, default_params, default_u0s, checkunits)
	structmodule = PlantModules.structmod(node)
	func_modules = [coupling.first for coupling in module_coupling if structmodule in coupling.second]

	component_systems = Vector{ODESystem}(undef, length(func_modules)) #! error if no func_modules found for structural module

	for (modulenum, func_module) in enumerate(func_modules)
		nodeparams = getnodeparamu0s(node, structmodule, func_module, module_defaults, default_params)
		nodeu0s = getnodeparamu0s(node, structmodule, func_module, module_defaults, default_u0s)

		component_systems[modulenum] = func_module(; :name => :foo,
			Pair.(keys(nodeparams), values(nodeparams))..., Pair.(keys(nodeu0s), values(nodeu0s))...)
				# real name given later
	end

	MTKsystem = collapse(
		component_systems,
		name = Symbol(string(structmodule) * string(PlantModules.id(node))),
		checkunits = checkunits
	)

	return MTKsystem
end

# get correct parameter/initial values for node between those defined in the model defaults, module defaults and node values
function getnodeparamu0s(node, structmodule, func_module, module_defaults, default_paramu0s)
	if !haskey(default_paramu0s, Symbol(func_module))
		return Dict() # no parameters specified in default_paramu0s => functional module has no parameters/initial values
	end

	paramu0s = default_paramu0s[Symbol(func_module)] |> x -> Dict{Any, Any}(Pair.(keys(x), values(x))) # PlantModules defaults
		#! excluding {Any, Any} causes the Dict to overspecialize on a type sometimes, better fix available?
	paramu0names = keys(paramu0s)

	node_attributes = PlantModules.attributes(node)
	nodemoduledefaults = module_defaults[structmodule]

	for paramu0name in paramu0names
		if paramu0name in keys(node_attributes) # Change to node-specific value
			paramu0s[paramu0name] = node_attributes[paramu0name]
		elseif paramu0name in keys(nodemoduledefaults) # Change to module-wide defaults
			paramu0s[paramu0name] = nodemoduledefaults[paramu0name]
		end
	end

	return paramu0s
end

# collapse multiple ODESystems into one. like ModelingToolkit.compose, but keeps a single namespace
function collapse(systems::Vector{ODESystem}; name::Symbol, checkunits::Bool)
    return ODESystem(
		vcat([get_eqs(system) for system in systems]...),
		get_iv(systems[1]), vcat([unknowns(system) for system in systems]...),
        vcat([parameters(system) for system in systems]...),
		name = name, checks = checkunits
		)
end #! use `extend` instead?

# Get the MTK system of the edge between the two nodes, and whether it exists in correct order
function get_connecting_module(node, nb_node, connecting_modules)
	structmodule = PlantModules.structmod(node)
	nb_structmodule = PlantModules.structmod(nb_node)

	connecting_module_idx = findfirst(x -> issetequal(x.first, (structmodule, nb_structmodule)), connecting_modules)

	if isnothing(connecting_module_idx)
		error("No connection module found for edges between nodes of types $structmodule and $nb_structmodule.")
	end

	connecting_module = connecting_modules[connecting_module_idx]
	reverse_order = connecting_module.first[1] == nb_structmodule
	return connecting_module, reverse_order
end

# get MTK system of connection between a node and its neighbour node AND the equations connecting the edge with the nodes
function get_connection_info(node, graphnr, nb_node, nb_node_graphnr, connecting_module, reverse_order,
	default_params, default_u0s, MTK_system_dicts
	)

	structmodule = PlantModules.structmod(node)
	nb_structmodule = PlantModules.structmod(nb_node)

	connector_func, connection_specific_values = connecting_module.second
	
	default_conn_info = merge( # merge parameter and initial values. #! if the default values are empty, assume there are none?
		get(default_params, Symbol(connector_func), ()),
		get(default_u0s, Symbol(connector_func), ())
	)
	conn_info = merge(default_conn_info, connection_specific_values)

	connection_MTK, get_connection_eqset = connector_func(;
		name = Symbol(string(structmodule) * string(PlantModules.id(node)) * "_" *
			string(nb_structmodule) * string(PlantModules.id(nb_node))),
		Pair.(keys(conn_info), values(conn_info))...
	)

	node_MTK = MTK_system_dicts[graphnr][PlantModules.id(node)]
	nb_node_MTK = MTK_system_dicts[nb_node_graphnr][PlantModules.id(nb_node)]

	connection_eqset = get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order)

	return connection_MTK, connection_eqset
end

# get neighbouring nodes of a node both from the same graph and all connected graphs
function get_nb_nodes(node, graphs, graphnr, intergraph_connections)
	graph = graphs[graphnr]
	
	intra_nb_nodes = PlantModules.neighbours(node, graph)
	intra_nb_node_graphnrs = repeat([graphnr], length(intra_nb_nodes))

	inter_nb_nodes, inter_nb_node_graphnrs = get_intergraph_neighbours(node, graphnr, graphs, intergraph_connections)

	nb_nodes = vcat(intra_nb_nodes, inter_nb_nodes)
	nb_node_graphnrs = vcat(intra_nb_node_graphnrs, inter_nb_node_graphnrs)

	return nb_nodes, nb_node_graphnrs
end

# get neighbouring nodes from different graphs
## go over all intergraph connections
function get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections::Vector)

	# Get intergraph connections involving the graph `node` belongs to
	node_intergraph_connections_plus_igidxs = [(intergraph_connection, findfirst(node_graphnr .== intergraph_connection[1]))
		for intergraph_connection in intergraph_connections 
		if node_graphnr in intergraph_connection[1]
	] # igidx is index of the node's graph in the intergraph connection

	nb_nodes = []
	nb_node_graphnrs = []
	
	# Go over all connected graphs and add nodes to neighbouring nodes if the specified condition is met
	for (node_intergraph_connection, igidx) in node_intergraph_connections_plus_igidxs
		nb_igidx = 3-igidx # 2 if igidx == 1 and 1 if igidx == 2
		nb_graphnr = node_intergraph_connection[1][nb_igidx]
		nb_graph = graphs[nb_graphnr]
		connection = node_intergraph_connection[2]
		if connection isa Tuple
			_nb_nodes = _get_intergraph_neighbours(node, nb_graph, connection[igidx], connection[nb_igidx])
		else
			nb_first = nb_igidx == 1
			_nb_nodes = _get_intergraph_neighbours(node, nb_graph, connection, nb_first)
		end
		
		append!(nb_nodes, _nb_nodes)
		append!(nb_node_graphnrs, repeat([nb_graphnr], length(_nb_nodes)))
	end

	return nb_nodes, nb_node_graphnrs
end

## For a connection where the structural module or specific nodes are specified
function _get_intergraph_neighbours(node, nb_graph, node_connection, nb_connection)
	if connection_check(node, node_connection)
		nb_nodes = [nb_node for nb_node in PlantModules.nodes(nb_graph) if connection_check(nb_node, nb_connection)]
	else
		nb_nodes = []
	end

	return nb_nodes
end

connection_check(node, connection) = node == connection # connection is a node
connection_check(node, connection::Symbol) = PlantModules.structmod(node) == connection # connection is a structural module
connection_check(node, connection::Vector) = node in connection # connection is a collection of nodes

## For a connection with a user-defined filter function
function _get_intergraph_neighbours(node, nb_graph, connection_func::Function, nb_first::Bool)
	if nb_first
		nb_nodes = [nb_node for nb_node in PlantModules.nodes(nb_graph) if connection_func(nb_node, node)]
	else
		nb_nodes = [nb_node for nb_node in PlantModules.nodes(nb_graph) if connection_func(node, nb_node)]
	end

	return nb_nodes
end