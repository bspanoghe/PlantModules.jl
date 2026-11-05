"""
    remake_graphsystem(prob::AbstractSciMLProblem, sys::System, structure, varnames, subsystem_types, value)

Remake the given `prob`, changing the values of subsystem variables to a specified value.
Only variables of given names and subsystem types are changed.

# Inputs
- `prob::AbstractSciMLProblem`: A SciML problem.
- `sys::System`: The ModelingToolkit.jl system corresponding to the given problem.
- `structure`: A graph representing the subsystem structure of the system. PlantModules' graph functions must be extended for the graph type. See also [`PlantStructure`](@ref).
- `varnames`: The name(s) of the desired variable.
- `subsystem_types`: The desired type(s) of subsystem. For node modules, this corresponds to a node or a structural module type. For edge modules (or connection modules), this corresponds to a connection, being a 2-tuple of node(s) and structural module(s).
- `value`: The new variable value.
"""
function remake_graphsystem(prob::AbstractSciMLProblem, sys::System, structure, varnames, subsystem_types, value)
    # can get even more efficient: https://docs.sciml.ai/ModelingToolkit/dev/examples/remake/
    remakevars = get_subsystem_variables(sys, structure, varnames, subsystem_types)
    ps = copy(parameter_values(prob))
    setter = setp(prob, remakevars)
    setter(ps, fill(value, length(remakevars)))
    newprob = remake(prob, p = ps)
    return newprob
end

"""
    remake_graphsystem!(prob::AbstractSciMLProblem, sys::System, structure, varnames::Symbol, subsystem_types, value)

Mutating version of [`remake_graphsystem`](@ref).
"""
function remake_graphsystem!(prob::AbstractSciMLProblem, sys::System, structure, varnames::Symbol, subsystem_types, value)
    remakevars = get_subsystem_variables(sys, structure, varnames, subsystem_types)
    ps = parameter_values(prob)
    setter = setp(prob, remakevars)
    setter(ps, fill(value, length(remakevars)))
    return remake(prob, p = ps)
end