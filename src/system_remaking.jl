"""
    remake_graphsystem(prob::AbstractSciMLProblem, value, sys::System, structure, varname::Symbol, subsystem_type)

Remake the given `prob`, changing the values of subsystem variables to a specified value.
Only variables of a given name and subsystem type are changed.

# Inputs
- `prob::AbstractSciMLProblem`: A SciML problem.
- `sys::System`: The ModelingToolkit.jl system corresponding to the given problem.
- `structure`: A graph representing the subsystem structure of the system. PlantModules' graph functions must be extended for the graph type. See also [`PlantStructure`](@ref).
- `varname`: The name of the desired variable.
- `subsystem_type`: The desired type of subsystem. For node modules, this corresponds to a node or structural module. For edge modules (or connection modules), this corresponds to a connection, being a 2-tuple of node(s) and structural module(s).
- `value`: The new variable value.
"""
function remake_graphsystem(prob::AbstractSciMLProblem, sys::System, structure, varname::Symbol, subsystem_type, value)
        # can get even more efficient: https://docs.sciml.ai/ModelingToolkit/dev/examples/remake/
    remakevars = get_subsystem_variables(sys, structure, varname, subsystem_type)
    ps = copy(parameter_values(prob))
    setter = setp(prob, remakevars)
    setter(ps, fill(value, length(remakevars)))
    newprob = remake(prob, p = ps)
    return newprob
end

"""
remake_graphsystem!(prob::AbstractSciMLProblem, sys::System, structure, varname::Symbol, subsystem_type, value)

Mutating version of [`remake_graphsystem`](@ref).
"""
function remake_graphsystem!(prob::AbstractSciMLProblem, sys::System, structure, varname::Symbol, subsystem_type, value)
    remakevars = get_subsystem_variables(sys, structure, varname, subsystem_type)
    ps = parameter_values(prob)
    setter = setp(prob, remakevars)
    setter(ps, fill(value, length(remakevars)))
    remake(prob, p = ps)
end

"""
    get_subsystem_variables(sys::System, structure, varname::Symbol, subsystem_type)

Get the Symbolics representation of all variables (unknowns or parameters) of a system with a given subsystem structure, filtered by variable name and type of subsystem.
See [`remake_graphsystem`](@ref) for more information about the inputs.
"""
function get_subsystem_variables(sys::System, structure, varname::Symbol, subsystem_type)
    nodes = getnodes(structure)
    sysnames = getsysnames(nodes, subsystem_type, structure)

    subsystems = [getsubsystem(sys, sysname) for sysname in sysnames]
    subsys_vars = [getproperty(subsys, varname) for subsys in subsystems]
    return subsys_vars
end

function getsysnames(nodes, node, _)
    @assert (node in nodes)
    sysnames = [string(getstructmod(node)) * string(getid(node))]
    return sysnames
end

function getsysnames(nodes, structmod::Symbol, _)
    node_structmods = getstructmod.(nodes)
    is_valid_node = [node_structmod == structmod for node_structmod in node_structmods]
    if !any(is_valid_node)
        error("None of the structural module \"$(structmod)\" were found in the graph.")
    end
    sysnames = [string(getstructmod(node)) * string(getid(node)) for node in nodes[is_valid_node]]
    return sysnames
end

function getsysnames(nodes, connection::Tuple, structure)
    is_valid_node = [
        [connection_check(node, connection[i]) for node in nodes] # see `graph_types.jl` for `connection_check`
        for i in eachindex(connection)
    ]
    if !any(Iterators.flatten(is_valid_node))
        error("No nodes found in the graph that correspond to connection $connection.")
    end

    sysnames = [
        string(getstructmod(node1)) * string(getid(node1)) * "_" * 
            string(getstructmod(node2)) * string(getid(node2))
        for node1 in nodes[is_valid_node[1]]
        for node2 in nodes[is_valid_node[2]]
        if node2 in getneighbors(node1, structure)
    ]
    isempty(sysnames) && error("No nodes found in the graph that correspond to connection $connection.")
    return sysnames
end

function getsubsystem(sys::System, sysname)
    parentsystem = get_parent(sys) # system before simplification
    subsystem = [subsys for subsys in get_systems(parentsystem) if get_name(subsys) == Symbol(sysname)][1]
    
    return subsystem
end