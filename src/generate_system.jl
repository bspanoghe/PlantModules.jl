# # Types

"""
	PlantStructure

A container for the plant structure. Used in [`generate_system`](@ref).

# Fields
- `graphs::Vector`: A vector of graphs representing the plants and their environment.
- `intergraph_connections::Vector`: A vector of pairs between 2 indexes of above graphs and how they are linked.
"""
struct PlantStructure
	graphs::Vector
	intergraph_connections::Vector
end

"""
	PlantCoupling

A container for the coupling of plant structure and function. Used in [`generate_system`](@ref).

# Fields
- `module_coupling::Dict`: Coupling between functional and structural modules.
- `connecting_modules::Dict{Tuple{Symbol, Symbol}, Function}`: The `System` to use for edges between specified nodes.
- `connecting_eqs::Function`: A function returning a vector of equations linking a node and all its edges.
"""
struct PlantCoupling
	module_coupling::Dict
	connecting_modules::Dict{Tuple{Symbol, Symbol}, Function}
	connecting_eqs::Function
end

"""
	PlantParameters

A container for functional parameters. Used in [`generate_system`](@ref).

# Fields
- `default_values::Dict{Symbol, <:Any}`: The model-wide default values of parameters and initial values.
- `module_defaults::Dict{Symbol, Dict}`: Module-specific default values of parameters and initial values.
- `connection_values::Dict{Tuple{Symbol, Symbol}, Dict}`: Connection-specific values of parameters and initial values.
"""
struct PlantParameters
	default_values::Dict{Symbol, <:Any}
	module_defaults::Dict{Symbol, Dict}
	connection_values::Dict{Tuple{Symbol, Symbol}, Dict}
end

# # Constructors

"""
	PlantStructure(graph)

Constructor function for `PlantStructure` variables.
"""
function PlantStructure(graph)
	return PlantStructure([graph], [])
end

function PlantStructure(::AbstractArray)
	error("`intergraph_connections` must be specified when the structure consists of multiple graphs.")
end

"""
	PlantCoupling(; module_coupling::Dict, connecting_modules::Dict,
		connecting_eqs::Function = PlantModules.multi_connection_eqs)

Constructor function for `PlantCoupling` variables.
"""
function PlantCoupling(; module_coupling::Dict, connecting_modules::Dict,
	connecting_eqs::Function = PlantModules.multi_connection_eqs)
	
	return PlantCoupling(module_coupling, connecting_modules, connecting_eqs)
end

"""
	PlantParameters(; default_values = PlantModules.default_values, 
		module_defaults = Dict(), extra_defaults = Dict())

Constructor function for `PlantParameters` variables.

`default_changes` is a dictionary with values to add to `default_values`. Refer to the `PlantParameters` type for a description of the other inputs.
"""
function PlantParameters(; default_values::Dict = PlantModules.default_values, module_defaults::Dict = Dict(), 
	connection_values::Dict = Dict(), default_changes::Dict = typeof(default_values)())

	changed_defaults = merge(default_values, default_changes)
	
	return PlantParameters(changed_defaults, module_defaults, connection_values)
end

# # Functions

"""
	generate_system(plantstructure::PlantStructure, plantparams::PlantParameters,
		plantcoupling::PlantCoupling; checkunits::Bool = true)

Create a ModelingToolkit.jl [`System`](@ref) that describes the functional behaviour of the input graph structure.

# Arguments
- `plantstructure`: One or more graphs specifying how the structural modules are connected, alongside the connections between graphs. See [`PlantStructure`](@ref).
- `plantparams`: Additional functional information about the connections between structural modules. See [`PlantParameters`](@ref).
- `plantcoupling`: Coupling between functional and structural modules. See [`PlantCoupling`](@ref).
- `checkunits`: Should the model check the units of the given equations? Defaults to `true`.
"""
function generate_system(plantstructure::PlantStructure, plantcoupling::PlantCoupling, plantparams::PlantParameters;
		checkunits::Bool = true)

	MTK_system_dicts = get_MTK_system_dicts(plantstructure, plantparams, plantcoupling, checkunits) 
		# Vector of a dict per graph with (node id => node MTK system)
	MTK_systems = vcat(collect.(values.(MTK_system_dicts))...) # node MTK systems

	connection_MTKs = System[] # MTK systems of edges
	connection_eqsets = Equation[] # equations linking nodes with edges

	for (graphnr, graph) in enumerate(plantstructure.graphs) # go over all graphs
		for node in PlantModules.getnodes(graph) # go over every node in the graph
			
			node_id = PlantModules.getid(node)
			nb_nodes, nb_node_graphnrs = get_nb_nodes(node, graphnr, plantstructure) # collect neighbour nodes
			current_connection_MTKs = Vector{System}(undef, length(nb_nodes)) # MTKs of current node's connections

			for (nb_idx, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs)) # go over all neighbours of the node
				connecting_module, reverse_order = get_connecting_module(node, nb_node, plantcoupling)
				connection_MTK, connection_eqset = get_connection_info( 
					node, graphnr, nb_node, nb_node_graphnr, connecting_module,
					reverse_order, plantparams, MTK_system_dicts
				)
				current_connection_MTKs[nb_idx] = connection_MTK
				append!(connection_eqsets, connection_eqset)
			end
			append!(connection_MTKs, current_connection_MTKs)
			append!(connection_eqsets, plantcoupling.connecting_eqs(MTK_system_dicts[graphnr][node_id], current_connection_MTKs))
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
function get_MTK_system_dicts(plantstructure, plantparams, plantcoupling, checkunits) 
	return [
		[PlantModules.getid(node) => 
			getMTKsystem(node, plantparams, plantcoupling, checkunits)
			for node in PlantModules.getnodes(graph)
		] |> Dict for graph in plantstructure.graphs
	]
end

# Get MTK system corresponding with node
function getMTKsystem(node, plantparams, plantcoupling, checkunits)
	structmodule = PlantModules.getstructmod(node)
	(!haskey(plantcoupling.module_coupling, structmodule) || isempty(plantcoupling.module_coupling[structmodule])) &&
		error("No functional module defined for structural module `$structmodule`.")
	func_modules = plantcoupling.module_coupling[structmodule]

	component_systems = Vector{System}(undef, length(func_modules))

	for (modulenum, func_module) in enumerate(func_modules)
		nodevalues = getnodevalues(node, structmodule, func_module, plantparams)

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
function getnodevalues(node, structmodule, func_module, plantparams)

	node_defaults = get_func_defaults(plantparams.default_values, func_module)
	node_module_defaults = get(plantparams.module_defaults, structmodule, Dict())
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
function get_connecting_module(node, nb_node, plantcoupling)
	structmodule = PlantModules.getstructmod(node)
	nb_structmodule = PlantModules.getstructmod(nb_node)

	if haskey(plantcoupling.connecting_modules, (structmodule, nb_structmodule))
		connecting_module = plantcoupling.connecting_modules[(structmodule, nb_structmodule)]
		reverse_order = false
	elseif haskey(plantcoupling.connecting_modules, (nb_structmodule, structmodule))
		connecting_module = plantcoupling.connecting_modules[(nb_structmodule, structmodule)]
		reverse_order = true
	else
		error("No connection module found for edges between nodes of types $structmodule and $nb_structmodule.")
	end

	return connecting_module, reverse_order
end

# get MTK system of connection between a node and its neighbour node AND the equations connecting the edge with the nodes
function get_connection_info(node, graphnr, nb_node, nb_node_graphnr, connecting_module,
		reverse_order, plantparams, MTK_system_dicts)

	structmodule = PlantModules.getstructmod(node)
	nb_structmodule = PlantModules.getstructmod(nb_node)

	connection_specific_values = reverse_order ?
		get(plantparams.connection_values, (nb_structmodule, structmodule), Dict()) : 
		get(plantparams.connection_values, (structmodule, nb_structmodule), Dict())
	
	default_values_conn = get_func_defaults(plantparams.default_values, connecting_module)
	conn_info = merge(default_values_conn, connection_specific_values)

	connection_MTK, get_connection_eqset = connecting_module(;
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
function get_nb_nodes(node, graphnr, plantstructure)
	graph = plantstructure.graphs[graphnr]
	
	intra_nb_nodes = PlantModules.getneighbours(node, graph)
	intra_nb_node_graphnrs = repeat([graphnr], length(intra_nb_nodes))

	inter_nb_nodes, inter_nb_node_graphnrs = get_intergraph_neighbours(node, graphnr, plantstructure)

	nb_nodes = vcat(intra_nb_nodes, inter_nb_nodes)
	nb_node_graphnrs = vcat(intra_nb_node_graphnrs, inter_nb_node_graphnrs)

	isempty(nb_nodes) && error("No neighbours found for node $node.")

	return nb_nodes, nb_node_graphnrs
end

# get neighbouring nodes from different graphs
## go over all intergraph connections
function get_intergraph_neighbours(node, node_graphnr, plantstructure)
	nb_nodes = []
	nb_node_graphnrs = []
	
	# Go over all graphs and add nodes to neighbouring nodes if the graphs are connected
	for intergraph_connection in plantstructure.intergraph_connections
		if node_graphnr in first(intergraph_connection)
			node_idx, nb_idx = first(intergraph_connection)[1] == node_graphnr ? (1, 2) : (2, 1)
			nb_graphnr = first(intergraph_connection)[nb_idx]
			nb_graph = plantstructure.graphs[nb_graphnr]
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