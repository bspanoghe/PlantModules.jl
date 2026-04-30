"""
    get_subsystem_variables(sys::System, structure, varname::Symbol, subsystem_type)

Get the Symbolics representation of all variables (unknowns or parameters) of a system with a given subsystem structure, filtered by variable name and type of subsystem.
See [`remake_graphsystem`](@ref) for more information about the inputs.
"""
function get_subsystem_variables(sys::System, structure, varnames::Vector{Symbol}, subsystem_types::Vector)
    nodes = getnodes(structure)
    sysnames = [getsysnames(nodes, subsystem_type, structure) for subsystem_type in subsystem_types] |>
        x -> reduce(vcat, x)

    subsystems = [getsubsystem(sys, sysname) for sysname in sysnames]
    subsys_vars = [getproperty(subsys, varname) for subsys in subsystems for varname in varnames]
    return subsys_vars
end

get_subsystem_variables(sys, structure, varnames::Vector, subsystem_type) = get_subsystem_variables(sys, structure, varnames, [subsystem_type])
get_subsystem_variables(sys, structure, varname, subsystem_types::Vector) = get_subsystem_variables(sys, structure, [varname], subsystem_types)
get_subsystem_variables(sys, structure, varname, subsystem_type) = get_subsystem_variables(sys, structure, [varname], [subsystem_type])

getsysname(node) = string(getstructmod(node)) * string(getid(node))

function getsysnames(nodes, node, _)
    (node in nodes) || error("Node $node not found.")
    sysnames = [getsysname(node)]
    return sysnames
end

function getsysnames(nodes, structmod::Symbol, _)
    node_structmods = getstructmod.(nodes)
    is_valid_node = [node_structmod == structmod for node_structmod in node_structmods]
    if !any(is_valid_node)
        error("None of the structural module \"$(structmod)\" were found in the graph.")
    end
    sysnames = [getsysname(node) for node in nodes[is_valid_node]]
    return sysnames
end

function getsysnames(nodes, connection::Tuple, structure)
    is_valid_node = [
        [connection_check(node, connection[i]) for node in nodes] # see `plantstructure.jl` for `connection_check`
            for i in eachindex(connection)
    ]
    if !any(Iterators.flatten(is_valid_node))
        error("No nodes found in the graph that correspond to connection $connection.")
    end

    sysnames = [
        getsysname(node1) * "_" * getsysname(node2)
        for node1 in nodes[is_valid_node[1]]
        for node2 in nodes[is_valid_node[2]]
        if node2 in getneighbors(node1, structure)
    ]
    isempty(sysnames) && error("No nodes found in the graph that correspond to connection $connection.")
    return sysnames
end

function getsubsystem(sys::System, sysname)
    parentsystem = get_parent(sys) # system before simplification
    if !isnothing(get_parent(parentsystem))
        parentsystem = get_parent(parentsystem) #! you have to do this twice since MTKv11
    end
    subsystems = get_systems(parentsystem)
    subsystem = subsystems[findfirst(subsys -> get_name(subsys) == Symbol(sysname), subsystems)] 

    return subsystem
end
