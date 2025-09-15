# # Types & constructors

"""
	PlantFunctionality

A container for functional parameters to be used in the function `generate_system`.

# Fields
- `graphs`: A vector of graphs representing the plants and their environment.
- `intergraph_connections`: A vector of pairs between 2 indexes of above graphs and how they are linked.
"""
struct PlantStructure
	graphs::Vector
	intergraph_connections::Vector
end

PlantStructure(graph) = PlantStructure([graph], [])

"""
	PlantFunctionality

A container for functional parameters to be used in the function `generate_system`.

# Fields
- `default_values`: Model-wide default values of parameters and initial values.
- `module_defaults`: Module-specific default values of parameters and initial values.
- `connecting_modules`: The `System` to use for edges between specified nodes.
- `connecting_eqs`: A function returning a vector of equations linking a node and all its edges.
"""
struct PlantFunctionality
	default_values::Dict{Symbol, <:Any}
	module_defaults::Dict{Symbol, Dict}
	connecting_modules::Vector
	connecting_eqs::Function
end

"""
	PlantFunctionality(; default_values = PlantModules.default_values, module_defaults = Dict(),
	connecting_modules, connecting_eqs = PlantModules.multi_connection_eqs, extra_defaults = Dict()
	)

Constructor function for `PlantFunctionality` variables.

`default_changes` is a dictionary with values to add to `default_values`. Refer to the `PlantFunctionality` type for a description of the other inputs.
"""
function PlantFunctionality(; default_values::Dict = PlantModules.default_values, module_defaults::Dict = Dict(),
	connecting_modules::Vector,	connecting_eqs::Function = PlantModules.multi_connection_eqs, 
	default_changes::Dict = typeof(default_values)()
	)

	changed_defaults = merge(default_values, default_changes)
	
	return PlantFunctionality(changed_defaults, module_defaults, connecting_modules, connecting_eqs)
end

# # Functions

"""
	generate_system(default_params::NamedTuple, default_u0s::NamedTuple, module_defaults::NamedTuple,
		module_coupling::Vector, struct_connections::Vector, func_connections::Vector; checkunits::Bool)

Creates a ModelingToolkit.jl [`System`](@ref) that describes the functional behaviour of the input graph structure.

# Arguments
- `default_params`: Model-wide default parameter values.
- `default_u0s`: Model-wide default initial values.
- `module_defaults`: Module-specific default values of both parameters and initial values.
- `module_coupling`: Coupling between functional and structural modules.
- `struct_connections`: One or more graphs specifying how the structural modules are connected, alongside the connections between graphs.
- `func_connections`: Additional functional information about the connections between structural modules.
- `checkunits`: Should the model check the units of the given equations? Defaults to `true`.
"""
function generate_system(struct_connections::PlantStructure, func_connections::PlantFunctionality,
	module_coupling::Dict; checkunits::Bool = true
	)

	MTK_system_dicts = get_MTK_system_dicts(
		struct_connections.graphs, func_connections.default_values, func_connections.module_defaults, module_coupling, checkunits
	) # Vector of a dict per graph with (node id => node MTK system)
	MTK_systems = vcat(collect.(values.(MTK_system_dicts))...) # node MTK systems

	connection_MTKs = System[] # MTK systems of edges
	connection_eqsets = Equation[] # equations linking nodes with edges

	for (graphnr, graph) in enumerate(struct_connections.graphs) # go over all graphs
		for node in PlantModules.getnodes(graph) # go over every node in the graph
			
			node_id = PlantModules.getid(node)
			nb_nodes, nb_node_graphnrs = get_nb_nodes(node, struct_connections.graphs, graphnr, struct_connections.intergraph_connections) # collect neighbour nodes
			current_connection_MTKs = Vector{System}(undef, length(nb_nodes)) # MTKs of current node's connections

			for (nb_idx, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs)) # go over all neighbours of the node
				connecting_module, reverse_order = get_connecting_module(node, nb_node, func_connections.connecting_modules)
				connection_MTK, connection_eqset = get_connection_info( 
					node, graphnr, nb_node, nb_node_graphnr, connecting_module,
					reverse_order, func_connections.default_values, MTK_system_dicts
				)
				current_connection_MTKs[nb_idx] = connection_MTK
				append!(connection_eqsets, connection_eqset)
			end
			append!(connection_MTKs, current_connection_MTKs)
			append!(connection_eqsets, func_connections.connecting_eqs(MTK_system_dicts[graphnr][node_id], current_connection_MTKs))
		end
	end

	model = compose(
		System(connection_eqsets, get_iv(MTK_systems[1]), name = :system, checks = checkunits),
		MTK_systems..., connection_MTKs...
	) # combine all subsystems together with equations linking nodes with edges

	system = mtkcompile(model)

	return system
end

# get node idx => node MTK system
function get_MTK_system_dicts(graphs, default_values, module_defaults, module_coupling, checkunits) 
	return [
		[PlantModules.getid(node) => 
			getMTKsystem(node, default_values, module_defaults, module_coupling, checkunits)
			for node in PlantModules.getnodes(graph)
		] |> Dict for graph in graphs
	]
end

# Get MTK system corresponding with node
function getMTKsystem(node, default_values, module_defaults, module_coupling, checkunits)
	structmodule = PlantModules.getstructmod(node)
	!haskey(module_coupling, structmodule) && error("No functional module defined for structural module `$structmodule`.")
	func_modules = module_coupling[structmodule]

	component_systems = Vector{System}(undef, length(func_modules)) #! error if no func_modules found for structural module

	for (modulenum, func_module) in enumerate(func_modules)
		nodevalues = getnodevalues(node, structmodule, func_module, module_defaults, default_values)

		component_systems[modulenum] = func_module(; :name => Symbol(string(structmodule) * string(PlantModules.getid(node))),
			Pair.(keys(nodevalues), values(nodevalues))...)
	end

	MTKsystem = component_systems[1]

	for comp_sys in component_systems[2:end]
		MTKsystem = extend(MTKsystem, comp_sys, checkunits)
	end

	return MTKsystem
end

# get correct parameter/initial values for node between those defined in the model defaults, module defaults and node values
function getnodevalues(node, structmodule, func_module, module_defaults, default_values)

	node_defaults = get_func_defaults(default_values, func_module)
	node_module_defaults = get(module_defaults, structmodule, Dict())
	node_attributes = PlantModules.getattributes(node)

	nodevalues = overwrite!(deepcopy(node_defaults), node_module_defaults, node_attributes)
	return nodevalues
end

# filters the values required for `funcmod` out of `default_values`
# intended to use with `System` functions so filters out `name` kwarg
function get_func_defaults(default_values::Dict, funcmod::Function)
	func_argnames = only(methods(funcmod)) |> Base.kwarg_decl |> x -> x[x .!= :name]
	return Dict([func_argname => default_values[func_argname] for func_argname in func_argnames if func_argname in keys(default_values)])
end

# overwrites values of first Dict with those in following Dicts 
function overwrite!(dicts::Dict...)
	commontype = promote_type(valtype.(dicts)...)
	maindict = convert(Dict{Symbol, commontype}, dicts[1])
	for dict in dicts[2:end]
		for key in keys(dict)
			if haskey(maindict, key)
				maindict[key] = dict[key]
			end
		end
	end

	return maindict
end

# extended version of ModelingToolkit.extend to include unitful checks
function extend(sys::ModelingToolkit.AbstractSystem, basesys::ModelingToolkit.AbstractSystem, checkunits::Bool; name::Symbol = nameof(sys),
	gui_metadata = get_gui_metadata(sys))

	T = SciMLBase.parameterless_type(basesys)
	ivs = independent_variables(basesys)
	if !(sys isa T)
		if length(ivs) == 0
			sys = convert_system(T, sys)
		elseif length(ivs) == 1
			sys = convert_system(T, sys, ivs[1])
		else
			throw("Extending multivariate systems is not supported")
		end
	end

	eqs = union(get_eqs(basesys), get_eqs(sys))
	sts = union(get_unknowns(basesys), get_unknowns(sys))
	ps = union(get_ps(basesys), get_ps(sys))
	base_deps = get_parameter_dependencies(basesys)
	deps = get_parameter_dependencies(sys)
	dep_ps = isnothing(base_deps) ? deps :
			isnothing(deps) ? base_deps : union(base_deps, deps)
	obs = union(get_observed(basesys), get_observed(sys))
	cevs = union(get_continuous_events(basesys), get_continuous_events(sys))
	devs = union(get_discrete_events(basesys), get_discrete_events(sys))
	defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
	syss = union(get_systems(basesys), get_systems(sys))

	if length(ivs) == 0
		T(eqs, sts, ps, observed = obs, defaults = defs, name = name, systems = syss,
			continuous_events = cevs, discrete_events = devs, gui_metadata = gui_metadata,
			parameter_dependencies = dep_ps, checks = checkunits)
	elseif length(ivs) == 1
		T(eqs, ivs[1], sts, ps, observed = obs, defaults = defs, name = name,
			systems = syss, continuous_events = cevs, discrete_events = devs,
			gui_metadata = gui_metadata, parameter_dependencies = dep_ps, 
			checks = checkunits)
	end
end

# Get the MTK system of the edge between the two nodes, and whether it exists in correct order
function get_connecting_module(node, nb_node, connecting_modules)
	structmodule = PlantModules.getstructmod(node)
	nb_structmodule = PlantModules.getstructmod(nb_node)

	connecting_module_idx = findfirst(x -> issetequal(x.first, (structmodule, nb_structmodule)), connecting_modules)

	if isnothing(connecting_module_idx)
		error("No connection module found for edges between nodes of types $structmodule and $nb_structmodule.")
	end

	connecting_module = connecting_modules[connecting_module_idx]
	reverse_order = connecting_module.first[1] == nb_structmodule
	return connecting_module, reverse_order
end

# get MTK system of connection between a node and its neighbour node AND the equations connecting the edge with the nodes
function get_connection_info(node, graphnr, nb_node, nb_node_graphnr, connecting_module,
	reverse_order, default_values, MTK_system_dicts
	)

	structmodule = PlantModules.getstructmod(node)
	nb_structmodule = PlantModules.getstructmod(nb_node)

	connector_func, connection_specific_values = connecting_module.second
	
	default_values_conn = get_func_defaults(default_values, connector_func)
	conn_info = merge(default_values_conn, connection_specific_values)

	connection_MTK, get_connection_eqset = connector_func(;
		name = Symbol(string(structmodule) * string(PlantModules.getid(node)) * "_" *
			string(nb_structmodule) * string(PlantModules.getid(nb_node))),
		Pair.(keys(conn_info), values(conn_info))...
	)

	node_MTK = MTK_system_dicts[graphnr][PlantModules.getid(node)]
	nb_node_MTK = MTK_system_dicts[nb_node_graphnr][PlantModules.getid(nb_node)]

	if applicable(get_connection_eqset, node_MTK, nb_node_MTK, connection_MTK) 
		# check if `reverse_order` is specified by user (needed for asymmetrical connections)
		connection_eqset = get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK)
	else
		connection_eqset = get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order)
	end

	return connection_MTK, connection_eqset
end

# get neighbouring nodes of a node both from the same graph and all connected graphs
function get_nb_nodes(node, graphs, graphnr, intergraph_connections)
	graph = graphs[graphnr]
	
	intra_nb_nodes = PlantModules.getneighbours(node, graph)
	intra_nb_node_graphnrs = repeat([graphnr], length(intra_nb_nodes))

	inter_nb_nodes, inter_nb_node_graphnrs = get_intergraph_neighbours(node, graphnr, graphs, intergraph_connections)

	nb_nodes = vcat(intra_nb_nodes, inter_nb_nodes)
	nb_node_graphnrs = vcat(intra_nb_node_graphnrs, inter_nb_node_graphnrs)

	#! add error message if no neighbours found

	return nb_nodes, nb_node_graphnrs
end

# get neighbouring nodes from different graphs
## go over all intergraph connections
function get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections::Vector)
	nb_nodes = []
	nb_node_graphnrs = []
	
	# Go over all graphs and add nodes to neighbouring nodes if the graphs are connected
	for intergraph_connection in intergraph_connections
		if node_graphnr in first(intergraph_connection)
			node_idx, nb_idx = first(intergraph_connection)[1] == node_graphnr ? (1, 2) : (2, 1)
			nb_graphnr = first(intergraph_connection)[nb_idx]
			nb_graph = graphs[nb_graphnr]
			connection = intergraph_connection[2]
			if connection isa Tuple
				_nb_nodes = _get_intergraph_neighbours(node, nb_graph, connection[node_idx], connection[nb_idx])
			else
				nb_first = nb_idx == 1
				_nb_nodes = _get_intergraph_neighbours(node, nb_graph, connection, nb_first)
			end
			
			append!(nb_nodes, _nb_nodes)
			append!(nb_node_graphnrs, repeat([nb_graphnr], length(_nb_nodes)))
		end
	end

	return nb_nodes, nb_node_graphnrs
end

## For a connection where the structural module or specific nodes are specified
function _get_intergraph_neighbours(node, nb_graph, node_connection, nb_connection)
	if connection_check(node, node_connection)
		nb_nodes = [nb_node for nb_node in PlantModules.getnodes(nb_graph) if connection_check(nb_node, nb_connection)]
	else
		nb_nodes = []
	end

	return nb_nodes
end

connection_check(node, connection) = node == connection # connection is a node
connection_check(node, connection::Symbol) = PlantModules.getstructmod(node) == connection # connection is a structural module
connection_check(node, connection::Vector) = node in connection # connection is a collection of nodes

## For a connection with a user-defined filter function
function _get_intergraph_neighbours(node, nb_graph, connection_func::Function, nb_first::Bool)
	if nb_first
		nb_nodes = [nb_node for nb_node in PlantModules.getnodes(nb_graph) if connection_func(nb_node, node)]
	else
		nb_nodes = [nb_node for nb_node in PlantModules.getnodes(nb_graph) if connection_func(node, nb_node)]
	end

	return nb_nodes
end