# # Types

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
	PlantCoupling(; module_coupling::Dict, connecting_modules::Dict,
		connecting_eqs::Function = PlantModules.multi_connection_eqs)

Constructor function for `PlantCoupling` variables.

Has the following default `connecting_eqs` function, which connects a node with all its edges:
```
multi_connection_eqs(node_MTK, connection_MTKs) = [
	node_MTK.ΣF ~ sum([connection_MTK.F for connection_MTK in connection_MTKs])
]
```
This function simply states that a node's net incoming water flux equals the sum of all connected water flows.
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
- `plantstructure`: A graph specifying how the structural modules are connected.
- `plantparams`: Additional functional information about the connections between structural modules. See [`PlantParameters`](@ref).
- `plantcoupling`: Coupling between functional and structural modules. See [`PlantCoupling`](@ref).
- `checkunits`: Should the model check the units of the given equations? Defaults to `true`.
"""
function generate_system(plantstructure, plantcoupling::PlantCoupling, plantparams::PlantParameters;
		checkunits::Bool = true)

	MTK_system_dict = get_MTK_system_dict(plantstructure, plantparams, plantcoupling, checkunits) 
		# Vector of a dict (node id => node MTK system)
	MTK_systems = collect(values(MTK_system_dict)) # node MTK systems

	connection_MTKs = System[] # MTK systems of edges
	connection_eqsets = Equation[] # equations linking nodes with edges

	for node in PlantModules.getnodes(plantstructure) # go over every node in the graph
		nb_nodes = getneighbors(node, plantstructure) # collect neighbour nodes
		current_connection_MTKs = Vector{System}(undef, length(nb_nodes)) # MTKs of current node's connections

		for (nb_idx, nb_node) in enumerate(nb_nodes) # go over all neighbours of the node
			connecting_module, original_order = get_connecting_module(node, nb_node, plantcoupling)
			connection_MTK, connection_eqset = get_connection_info( 
				node, nb_node, connecting_module, original_order, plantparams, MTK_system_dict
			)
			current_connection_MTKs[nb_idx] = connection_MTK
			append!(connection_eqsets, connection_eqset)
		end
		append!(connection_MTKs, current_connection_MTKs)
		append!(connection_eqsets, plantcoupling.connecting_eqs(MTK_system_dict[PlantModules.getid(node)], current_connection_MTKs))
	end

	model = compose(
		System(connection_eqsets, get_iv(MTK_systems[1]), name = :system, checks = checkunits),
		MTK_systems..., connection_MTKs...
	) # combine all subsystems together with equations linking nodes with edges

	system = mtkcompile(model)

	return system
end

# get node idx => node MTK system
function get_MTK_system_dict(plantstructure, plantparams, plantcoupling, checkunits) 
	return [
		PlantModules.getid(node) => getMTKsystem(node, plantparams, plantcoupling, checkunits)
		for node in PlantModules.getnodes(plantstructure)
		] |> Dict
end

# Get MTK system corresponding with node
function getMTKsystem(node, plantparams, plantcoupling, checkunits)
	structmodule = PlantModules.getstructmod(node)
	(!haskey(plantcoupling.module_coupling, structmodule) || isempty(plantcoupling.module_coupling[structmodule])) &&
		error("No functional module defined for structural module `$structmodule`.")
	func_modules = plantcoupling.module_coupling[structmodule]

<<<<<<< HEAD
	component_systems = Vector{System}(undef, length(func_modules))
=======
	component_systems = Vector{ODESystem}(undef, length(func_modules))
>>>>>>> c71c803cbf5e84f69a5edcf0ed7a38b1cb5d4b7b

	for (modulenum, func_module) in enumerate(func_modules)
		nodevalues = getnodevalues(node, structmodule, func_module, plantparams)

<<<<<<< HEAD
		component_systems[modulenum] = func_module(; :name => Symbol(string(structmodule) * string(PlantModules.getid(node))),
			Pair.(keys(nodevalues), values(nodevalues))...)
=======
		component_systems[modulenum] = func_module(; :name => Symbol(string(structmodule) * string(PlantModules.id(node))),
			Pair.(keys(nodeparams), values(nodeparams))..., Pair.(keys(nodeu0s), values(nodeu0s))...)
				# real name given later
>>>>>>> c71c803cbf5e84f69a5edcf0ed7a38b1cb5d4b7b
	end

	MTKsystem = component_systems[1]

	for comp_sys in component_systems[2:end]
		MTKsystem = extend(MTKsystem, comp_sys, checkunits)
	end

	return MTKsystem
end

# extended version of ModelingToolkit.extend to include unitful checks
function extend(sys::AbstractSystem, basesys::AbstractSystem, checkunits::Bool; name::Symbol = nameof(sys),
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

# get correct parameter/initial values for node between those defined in the model defaults, module defaults and node values
function getnodevalues(node, structmodule, func_module, plantparams)
	node_defaults = get_func_defaults(plantparams.default_values, func_module)
	node_module_defaults = get(plantparams.module_defaults, structmodule, Dict())
	node_attributes = PlantModules.getattributes(node)

<<<<<<< HEAD
	nodevalues = overwrite!(deepcopy(node_defaults), node_module_defaults, node_attributes)
	return nodevalues
end
=======
	paramu0s = default_paramu0s[Symbol(func_module)] |> x -> Dict{Any, Any}(Pair.(keys(x), values(x))) # PlantModules defaults
	paramu0names = keys(paramu0s)
>>>>>>> c71c803cbf5e84f69a5edcf0ed7a38b1cb5d4b7b

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

<<<<<<< HEAD
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

=======
>>>>>>> c71c803cbf5e84f69a5edcf0ed7a38b1cb5d4b7b
# Get the MTK system of the edge between the two nodes, and whether it exists in correct order
function get_connecting_module(node, nb_node, plantcoupling)
	structmodule = PlantModules.getstructmod(node)
	nb_structmodule = PlantModules.getstructmod(nb_node)

	if haskey(plantcoupling.connecting_modules, (structmodule, nb_structmodule))
		connecting_module = plantcoupling.connecting_modules[(structmodule, nb_structmodule)]
		original_order = true
	elseif haskey(plantcoupling.connecting_modules, (nb_structmodule, structmodule))
		connecting_module = plantcoupling.connecting_modules[(nb_structmodule, structmodule)]
		original_order = false
	else
		error("No connection module found for edges between nodes of types $structmodule and $nb_structmodule.")
	end

	return connecting_module, original_order
end

# get MTK system of connection between a node and its neighbour node AND the equations connecting the edge with the nodes
function get_connection_info(node, nb_node, connecting_module,
		original_order, plantparams, MTK_system_dict)

	structmodule = PlantModules.getstructmod(node)
	nb_structmodule = PlantModules.getstructmod(nb_node)

	connection_specific_values = original_order ?
		get(plantparams.connection_values, (structmodule, nb_structmodule), Dict()) :
		get(plantparams.connection_values, (nb_structmodule, structmodule), Dict())
	
<<<<<<< HEAD
	default_values_conn = get_func_defaults(plantparams.default_values, connecting_module)
	order_info = :original_order in (connecting_module |> methods |> only |> Base.kwarg_decl) ?
		Dict{Symbol, Any}(:original_order => original_order) : Dict()
=======
	default_conn_info = merge(
		get(default_params, Symbol(connector_func), ()),
		get(default_u0s, Symbol(connector_func), ())
	)
	conn_info = merge(default_conn_info, connection_specific_values)
>>>>>>> c71c803cbf5e84f69a5edcf0ed7a38b1cb5d4b7b

	conn_info = merge(default_values_conn, connection_specific_values, order_info)

	connection_MTK, get_connection_eqset = connecting_module(;
		name = Symbol(string(structmodule) * string(PlantModules.getid(node)) * "_" *
			string(nb_structmodule) * string(PlantModules.getid(nb_node))),
		Pair.(keys(conn_info), values(conn_info))...
	)

	node_MTK = MTK_system_dict[PlantModules.getid(node)]
	nb_node_MTK = MTK_system_dict[PlantModules.getid(nb_node)]

	if applicable(get_connection_eqset, node_MTK, nb_node_MTK, connection_MTK) 
		# check if `original_order` is specified by user (needed for asymmetrical connections)
		connection_eqset = get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK)
	else
		connection_eqset = get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, original_order)
	end

	return connection_MTK, connection_eqset
end