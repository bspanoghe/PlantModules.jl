#! Maybe add RecipesBase.jl (get rid of Plots dependency)?

# Plot MTK solutions #

## Get the symbol representation of a MTK unknown (variable)
get_MTKunknown_symbol(s::SymbolicUtils.Symbolic) = s.metadata[ModelingToolkit.VariableSource][2]

## Pry the ODE system corresponding with given node out of the ODE solution
function getnodesystem(sol::ODESolution, node)
    system = getfield(sol.prob.f, :sys)
    nodename = string(PlantModules.structmod(node)) * string(PlantModules.id(node))
    nodesystem = [subsys for subsys in getproperty(system, :systems) if get_name(subsys) == Symbol(nodename)][1]
    return nodesystem
end

## Get all ODE systems corresponding with nodes off a certain structural module from an ODE solution
function getnodesystems(sol::ODESolution, graphnodes, struct_module::Symbol)
    matching_nodes = [node for node in graphnodes if PlantModules.structmod(node) == struct_module]
    if isempty(matching_nodes)
        error("No structural module $(struct_module) found in graph")
    end
    return [getnodesystem(sol, node) for node in matching_nodes]
end

get_graphnodes(graph) = PlantModules.nodes(graph)
get_graphnodes(graph::Vector) = reduce(vcat, get_graphnodes.(graph))

"""
    plotnode(sol::ODESolution, node; func_varname::Symbol)

Returns a plot for every functional variable of the given node for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable.
"""
function plotnode(sol::ODESolution, node; func_varname::Symbol = Symbol(""), kwargs...)
    nodesystem = getnodesystem(sol, node)
    indep_values = copy(sol[independent_variable(sol.prob.f.sys)]) # values of indepedent variable

    if func_varname == Symbol("") # No functional variable name given => show all
        func_varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystem)] |> unique
    else
        func_varnames = [func_varname]
    end

    var_values_vec = [sol[getproperty(nodesystem, func_varname)] for func_varname in func_varnames]
    plots = Vector{Plots.Plot}(undef, length(func_varnames))

    for (i, (var_values, func_varname)) in enumerate(zip(var_values_vec, func_varnames))
        if var_values isa Vector # multidimensional variable (e.g. vector of dimensions D)
            plots[i] = plot(indep_values, reduce(hcat, var_values)', title = "$func_varname"; kwargs...)
        else
            plots[i] = plot(indep_values, var_values, title = "$func_varname"; kwargs...)
        end
    end
  
    return only(plots)
end

"""
    plotgraph(sol::ODESolution, graph; struct_module::Symbol, func_varname::Symbol)

Returns a plot for every functional variable for every node of the given graph for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable,
give the name of a structural module to limit considered nodes to those of this type, or both.
"""

function plotgraph(sol::ODESolution, graph; struct_module::Symbol = Symbol(""), func_varname::Symbol = Symbol(""), kwargs...)

    indep_values = copy(sol[independent_variable(sol.prob.f.sys)]) # values of indepedent variable
    plot_data = Dict{Symbol, Dict{Symbol, Vector{Float64}}}()
    graphnodes = get_graphnodes(graph)

    if struct_module == Symbol("")
        struct_modules = graphnodes .|> PlantModules.structmod |> unique
    else
        struct_modules = [struct_module]
    end

    for struct_module in struct_modules
        nodesystems = getnodesystems(sol, graphnodes, struct_module)

        if func_varname == Symbol("")
            func_varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[1])] |> unique
        else
            func_varnames = [func_varname]
        end

        for func_varname in func_varnames
            if !haskey(plot_data, func_varname)
                plot_data[func_varname] = Dict{Symbol, Vector{Float64}}()
            end

            var_data = get!(plot_data[func_varname], struct_module, Float64[])

            for nodesystem in nodesystems
                var_values = sol[getproperty(nodesystem, func_varname)]

                if var_values[1] isa Vector # multidimensional variable (e.g. vector of dimensions D)
                    for var_dimension in eachrow(reduce(hcat, var_values))
                        append!(var_data, var_dimension)
                    end
                else
                    append!(var_data, var_values)
                end

            end
        end
    end

    plots = []

    for func_varname in keys(plot_data)
        var_data = plot_data[func_varname] |> values |> collect
        group_sizes = length.(var_data)
        groups = [fill(structmod, group_size) for (structmod, group_size) in zip(keys(plot_data[func_varname]), group_sizes)] |> 
            x -> reduce(vcat, x)

        ys = reduce(vcat, var_data)
        xs = repeat(indep_values, length(ys) ÷ length(indep_values))
        push!(plots, plot(xs, ys, group = groups, title = "$func_varname"; kwargs...))
    end

    return only(plots)
end