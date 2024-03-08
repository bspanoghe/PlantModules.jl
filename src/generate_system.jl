"""
	generate_system(model_defaults, module_defaults,	module_coupling,
	struct_connections, func_connections; checkunits = true)

Creates a MTK system based on a set of structural and functional modules and how those are connected.

# Arguments
- `model_defaults`: Model-wide default parameters.
- `module_defaults`: Module-specific default parameters.
- `module_coupling`: Coupling between functional and structural modules.
- `struct_connections`: One or more graphs specifying how the structural modules are connected, alongside the connections between graphs.
- `func_connections`: Additional functional information about the connections between structural modules.
- `checkunits`: Should the model check the units of the given equations? Defaults to `true`.
"""
function generate_system(model_defaults, module_defaults, module_coupling,
	struct_connections, func_connections; default_params = PlantModules.default_params,
	default_u0s = PlantModules.default_u0s, checkunits = true
	)

	graphs, intergraph_connections = struct_connections
    connecting_modules, get_connection_eqs = func_connections

	MTK_system_dicts = get_MTK_system_dicts(graphs, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
	MTK_systems = vcat(collect.(values.(MTK_system_dicts))...)

	connection_MTKs = ODESystem[]
	connection_equations = Equation[]

	for (graphnr, graph) in enumerate(graphs)
		for node in PlantModules.nodes(graph)
			nb_nodes, nb_node_graphnrs = get_nb_nodes(node, graphs, graphnr, intergraph_connections)

			if !isempty(nb_nodes) # no neighbours no connection info
				node_connection_MTK, node_connection_equations = get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, connecting_modules, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
				append!(connection_MTKs, node_connection_MTK)
				append!(connection_equations, node_connection_equations)
			end
		end
	end

	system = ODESystem(connection_equations, independent_variable(MTK_systems[1]), name = :system, systems = vcat(MTK_systems, connection_MTKs), checks = checkunits)

	return system
end

# get node idx to corresponding MTK system
get_MTK_system_dicts(graphs, module_coupling, module_defaults, model_defaults, default_params, default_u0s) = 
	[
		[PlantModules.id(node) => 
			getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
			for node in PlantModules.nodes(graph)
		] |> Dict for graph in graphs
]

# Get MTK system corresponding with node
function getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
	structmodule = PlantModules.structmod(node)
	func_modules = [coupling.first for coupling in module_coupling if structmodule in coupling.second]

	component_systems = Vector{ODESystem}(undef, length(func_modules))

	for (modulenum, func_module) in enumerate(func_modules)
		nodeparams = getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params)
		nodeu0s = getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)

		component_systems[modulenum] = func_module(; :name => :foo,
			Pair.(keys(nodeparams), values(nodeparams))..., Pair.(keys(nodeu0s), values(nodeu0s))...)
				# real name given later
	end

	MTKsystem = collapse(component_systems, name = Symbol(string(structmodule) * string(PlantModules.id(node))))

	return MTKsystem
end

# get correct parameter values for node between those defined in the model defaults, module defaults and node values
function getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params)
	if !haskey(default_params, Symbol(func_module))
		return Dict() # no parameters specified in default_params => functional module has no parameters
	end

	params = default_params[Symbol(func_module)] |> x -> Dict{Any, Any}(Pair.(keys(x), values(x))) # PlantModules defaults
		#! excluding {Any, Any} causes the Dict to overspecialize on a type sometimes, better fix available?
	paramnames = keys(params)

	node_attributes = PlantModules.attributes(node)
	nodemoduledefaults = module_defaults[structmodule]

	for paramname in paramnames
		if paramname in keys(node_attributes) # Change to node-specific value
			params[paramname] = node_attributes[paramname]
		elseif paramname in keys(nodemoduledefaults) # Change to module-wide defaults
			params[paramname] = nodemoduledefaults[paramname]
		elseif paramname in keys(model_defaults) # Change to model-wide defaults
			params[paramname] = model_defaults[paramname]
		end
	end

	return params
end

# get correct inital variable values for node between those defined in the model defaults, module defaults and node values
function getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)
	if !haskey(default_u0s, Symbol(func_module))
		return Dict() # no u0s specified in default_u0s => functional module has no u0s
	end
	
	u0s = default_u0s[Symbol(func_module)] |> x -> Dict{Any, Any}(Pair.(keys(x), values(x))) # PlantModules defaults
		#! excluding {Any, Any} causes the Dict to overspecialize on a type sometimes, better fix available?
	u0names = keys(u0s)

	node_attributes = PlantModules.attributes(node)
	nodemoduledefaults = module_defaults[structmodule]

	for u0name in u0names
		if u0name in keys(node_attributes) # Change to node-specific value
			u0s[u0name] = node_attributes[u0name]
		elseif u0name in keys(nodemoduledefaults) # Change to module-wide defaults
			u0s[u0name] = nodemoduledefaults[u0name]
		elseif u0name in keys(model_defaults) # Change to model-wide defaults
			u0s[u0name] = model_defaults[u0name]
		end
	end

	return u0s
end

# collapse multiple ODESystems into one. like ModelingToolkit.compose, but keeps a single namespace
function collapse(systems::Vector{ODESystem}; name::Symbol)
    return ODESystem(vcat([system.eqs for system in systems]...), systems[1].iv, vcat([unknowns(system) for system in systems]...),
        vcat([parameters(system) for system in systems]...), name = name)
end

# get MTK systems and the vector of equations defining connections between MTK systems of a node and its neighbour nodes
function get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, connecting_modules, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
	node_MTK = MTK_system_dicts[graphnr][PlantModules.id(node)]

	nb_node_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
	connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes))

	for (nb_nr, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs))
		nb_node_MTKs[nb_nr] = MTK_system_dicts[nb_node_graphnr][PlantModules.id(nb_node)]
		connection_MTKs[nb_nr] = get_func_connection(node, nb_node, connecting_modules, default_params, default_u0s)
	end

	connection_equations = get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs)
	return connection_MTKs, connection_equations
end

# given two nodes' structural modules, get the MTK system of the functional connection between them
function get_func_connection(node, nb_node, connecting_modules, default_params, default_u0s)
	structmodule = PlantModules.structmod(node)
	nb_structmodule = PlantModules.structmod(nb_node)

	connecting_module_vec = [connecting_module for connecting_module in connecting_modules if issetequal(connecting_module.first, (structmodule, nb_structmodule))]

	if isempty(connecting_module_vec)
		default_connecting_module = [connecting_module for connecting_module in connecting_modules if isempty(connecting_module.first)]
		connector_func, connection_specific_values = default_connecting_module[1].second, []
	else
		connector_func, connection_specific_values = connecting_module_vec[1].second
	end

	default_conn_info = merge(get(default_params, Symbol(connector_func), []), get(default_u0s, Symbol(connector_func), []))
	conn_info = merge(default_conn_info, connection_specific_values)

	func_connection = connector_func(;
		name = Symbol(string(structmodule) * string(PlantModules.id(node)) * "_" *
			string(nb_structmodule) * string(PlantModules.id(nb_node))),
		Pair.(keys(conn_info), values(conn_info))...
	)

	return func_connection
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

function get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
	nb_nodes = []
	nb_node_graphnrs = []
	
	node_intergraph_connections_plus_igidxs = [(intergraph_connection, findfirst(node_graphnr .== intergraph_connection[1]))
		for intergraph_connection in intergraph_connections 
		if node_graphnr in intergraph_connection[1]
	] # igidx is index of the node's graph in the intergraph connection

	for (node_intergraph_connection, igidx) in node_intergraph_connections_plus_igidxs
		if PlantModules.structmod(node) == node_intergraph_connection[2][igidx]
			nb_igidx = 3-igidx # 2 if igidx == 1 and 1 if igidx == 2
			nb_graphnr = node_intergraph_connection[1][nb_igidx]
			nb_graph = graphs[nb_graphnr]
			connection_nb_nodes = [nb_node for nb_node in PlantModules.nodes(nb_graph) if PlantModules.structmod(nb_node) == node_intergraph_connection[2][nb_igidx]]
			append!(nb_nodes, connection_nb_nodes)
			append!(nb_node_graphnrs, repeat([nb_graphnr], length(connection_nb_nodes)))
		end
	end

	return nb_nodes, nb_node_graphnrs
end