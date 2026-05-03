"""
    get_subsystem_variables(sys::System, plantstructure::PlantStructure, varnames::Vector{Symbol}, nodes::Vector{<:PMVertex})

Get the Symbolics representation of the variables (unknowns or parameters) of a system for a set of nodes, filtered by variable name.
See [`remake_graphsystem`](@ref) for more information about the inputs.
"""
function get_subsystem_variables(sys::System, plantstructure::PlantStructure, varnames::Vector{Symbol}, nodes::Vector{<:PMVertex})
    sysnames = [getsysname(node, plantstructure) for node in nodes]
    return _get_subsystem_variables(sys, sysnames, varnames)
end

function get_subsystem_variables(sys::System, plantstructure::PlantStructure, varnames::Vector{Symbol}, structmods::Vector{Symbol})
    sysnames = [getsysname(node, plantstructure) for node in getnodes(plantstructure) if getstructmod(node) in structmods]
    return _get_subsystem_variables(sys, sysnames, varnames)
end

function get_subsystem_variables(sys::System, plantstructure::PlantStructure, varnames::Vector{Symbol}, connections::Vector{<:Tuple})
    sysnames = [
        getsysname(node1, plantstructure) * "_" * getsysname(node2, plantstructure)
        for node1 in getnodes(plantstructure) for node2 in getnodes(plantstructure)
        if node2 in getneighbors(node1, plantstructure) &&
            any([
                connection_check(node1, connection[1]) && connection_check(node2, connection[2]) 
                for connection in connections
            ])
    ] #! CHECK

    return _get_subsystem_variables(sys, sysnames, varnames)
end

get_subsystem_variables(sys::System, plantstructure::PlantStructure, varnames::Vector, subsystem_type) = get_subsystem_variables(sys, plantstructure, varnames, [subsystem_type])
get_subsystem_variables(sys::System, plantstructure::PlantStructure, varname, subsystem_types::Vector) = get_subsystem_variables(sys, plantstructure, [varname], subsystem_types)
get_subsystem_variables(sys::System, plantstructure::PlantStructure, varname, subsystem_type) = get_subsystem_variables(sys, plantstructure, [varname], [subsystem_type])

"""
    get_subsystem_variables(sys::System, structure, graph_idx = 1, varnames::Vector{Symbol}, nodes::Vector{<:PMVertex})

Alternate method for non-PlantStructure graphs. Also requires the index of the graph has in the vector `graphs` that is passed to `PlantStructure` (set to 1 if not applicable).
"""
function get_subsystem_variables(sys::System, structure, graph_idx::Integer, varnames::Vector{Symbol}, nodes::Vector)
    graph_idx == 0 && error("If input graph is not of type PlantStructure, give the index that the graph has in the vector `graphs` that is passed to `PlantStructure` as keyword argument `graph_idx` (set to 1 if it is the only graph).")
    sysnames = [getsysname(node, graph_idx) for node in nodes]
    return _get_subsystem_variables(sys, sysnames, varnames)
end

function get_subsystem_variables(sys::System, structure, graph_idx::Integer, varnames::Vector{Symbol}, structmods::Vector{Symbol})
    graph_idx == 0 && error("If input graph is not of type PlantStructure, give the index that the graph has in the vector `graphs` that is passed to `PlantStructure` as keyword argument `graph_idx` (set to 1 if it is the only graph).")
    sysnames = [getsysname(node, graph_idx) for node in getnodes(structure) if getstructmod(node) in structmods]
    return _get_subsystem_variables(sys, sysnames, varnames)
end

get_subsystem_variables(sys, structure, graph_idx::Integer, varnames::Vector{Symbol}, subsystem_type) = get_subsystem_variables(sys, structure, graph_idx, varnames, [subsystem_type])
get_subsystem_variables(sys, structure, graph_idx::Integer, varname::Symbol, subsystem_types::Vector) = get_subsystem_variables(sys, structure, graph_idx, [varname], subsystem_types)
get_subsystem_variables(sys, structure, graph_idx::Integer, varname::Symbol, subsystem_type) = get_subsystem_variables(sys, structure, graph_idx, [varname], [subsystem_type])




# internal version
function _get_subsystem_variables(sys::System, sysnames, varnames)
    subsystems = [getsubsystem(sys, sysname) for sysname in sysnames]
    subsys_vars = [getproperty(subsys, varname) for subsys in subsystems for varname in varnames]

    length(subsys_vars) == 1 && return only(subsys_vars)
    return subsys_vars
end

# get name of the MTK system for a given node
getsysname(node, plantstructure::PlantStructure) = (
    Symbol(string(getstructmod(node)) * string(og_id(plantstructure, node)[1]) * "_" * string(og_id(plantstructure, node)[2]))
)

getsysname(node, graph_idx::Integer) = (
    Symbol(string(getstructmod(node)) * string(graph_idx) * "_" * string(getid(node)))
)

# get names of the MTK systems for a collection of nodes (as defined as a vector of nodes, a structural module type or a connection type)

function getsysnames(nodes, structmod::Symbol, plantstructure::PlantStructure)
    node_structmods = getstructmod.(nodes)
    is_valid_node = [node_structmod == structmod for node_structmod in node_structmods]
    if !any(is_valid_node)
        error("None of the structural module \"$(structmod)\" were found in the graph.")
    end
    sysnames = [getsysname(node, plantstructure) for node in nodes[is_valid_node]]
    return sysnames
end

function getsysnames(nodes, connection::Tuple, plantstructure::PlantStructure)
    is_valid_node = [
        [connection_check(node, connection[i]) for node in nodes] # see `plantstructure.jl` for `connection_check`
            for i in eachindex(connection)
    ]
    if !any(Iterators.flatten(is_valid_node))
        error("No nodes found in the graph that correspond to connection $connection.")
    end

    sysnames = [
        getsysname(node1, plantstructure) * "_" * getsysname(node2, plantstructure)
        for node1 in nodes[is_valid_node[1]]
        for node2 in nodes[is_valid_node[2]]
        if node2 in getneighbors(node1, plantstructure)
    ]
    isempty(sysnames) && error("No nodes found in the graph that correspond to connection $connection.")
    return sysnames
end

function getsubsystem(sys::System, sysname)
    parentsystem = get_parent(sys) # system before simplification
    if !isnothing(get_parent(parentsystem))
        parentsystem = get_parent(parentsystem) # note: you have to do this twice since MTKv11
    end
    subsystems = get_systems(parentsystem)
    subsystem_idx = findfirst(subsys -> get_name(subsys) == sysname, subsystems)
    isnothing(subsystem_idx) && error("Subsystem $(sysname) not found.")
    subsystem = subsystems[subsystem_idx] 

    return subsystem
end

function hassubsystem(sys::System, sysname)
    parentsystem = get_parent(sys) # system before simplification
    if !isnothing(get_parent(parentsystem))
        parentsystem = get_parent(parentsystem) # note: you have to do this twice since MTKv11
    end
    subsystems = get_systems(parentsystem)
    has_subsystem = sysname in get_name.(subsystems)
    return has_subsystem
end