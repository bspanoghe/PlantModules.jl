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
function getnodesystems(sol::ODESolution, graphnodes, structmod::Symbol)
    matching_nodes = [node for node in graphnodes if PlantModules.structmod(node) == structmod]
    if isempty(matching_nodes)
        error("No structural module $(structmod) found in graph")
    end
    return [getnodesystem(sol, node) for node in matching_nodes]
end

get_graphnodes(graph) = PlantModules.nodes(graph)
get_graphnodes(graph::Vector) = reduce(vcat, get_graphnodes.(graph))

# get variable values from solution
function get_values(sol, nodesystem, varname)
    if (varname in [:M, :Π, :Ψ]) & (:M_value in Symbol.(nodesystem.ps))
        # if the system has the "M_value" parameter, variable M is constant
        # this means the system was not solved for M and indexing for it or variables depending on it is slow
        # circumvent this by calculating it here instead
        system_param_dict = nodesystem.defaults |> x -> [Symbol(key) => val for (key, val) in x] |> Dict
        M_value = system_param_dict[:M_value]
        if varname == :M
            return fill(M_value, length(sol.t))
        else
            R = 8.314
            Π = - R * system_param_dict[:T] * M_value
            if varname == :Π
                return fill(Π, length(sol.t))
            elseif varname == :Ψ
                P = sol[getproperty(nodesystem, :P)]
                return P .+ Π
            end
        end
    elseif (varname == :K) & (:K_value in Symbol.(nodesystem.ps))
        system_param_dict = nodesystem.defaults |> x -> [Symbol(key) => val for (key, val) in x] |> Dict
        K_value = system_param_dict[:K_value]
        return fill(K_value, length(sol.t))
    else
        return sol[getproperty(nodesystem, varname)]
    end
end

"""
    plotnode(sol::ODESolution, node; varname::Symbol)

Returns a plot for every functional variable of the given node for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable.
"""
function plotnode(sol::ODESolution, node; varname::Symbol = Symbol(""), kwargs...)
    nodesystem = getnodesystem(sol, node)
    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable

    if varname == Symbol("") # No functional variable name given => show all
        varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystem)] |> unique
    else
        varnames = [varname]
    end

    var_values_vec = [get_values(sol, nodesystem, varname) for varname in varnames]
    plots = Vector{Plots.Plot}(undef, length(varnames))

    for (i, (var_values, varname)) in enumerate(zip(var_values_vec, varnames))
        if var_values isa Vector # multidimensional variable (e.g. vector of dimensions D)
            plots[i] = plot(indep_values, reduce(hcat, var_values)', title = "$varname"; kwargs...)
        else
            plots[i] = plot(indep_values, var_values, title = "$varname"; kwargs...)
        end
    end
  
    return length(plots) == 1 ? only(plots) : plots
end

"""
    plotgraph(sol::ODESolution, graph; structmod::Symbol, varname::Symbol)

Returns a plot for every functional variable for every node of the given graph for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable,
give the name of a structural module to limit considered nodes to those of this type, or both.
"""
function plotgraph(sol::ODESolution, graph; structmod::Symbol = Symbol(""), varname::Symbol = Symbol(""), kwargs...)

    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
    append!(indep_values, NaN) # NaN used to cause linebreaks in plot
    plot_data = Dict{Symbol, Dict{Symbol, Vector{Float64}}}()
    graphnodes = get_graphnodes(graph)

    if structmod == Symbol("")
        structmods = graphnodes .|> PlantModules.structmod |> unique
    else
        structmods = [structmod]
    end

    for structmod in structmods
        nodesystems = getnodesystems(sol, graphnodes, structmod)

        if varname == Symbol("")
            varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[1])] |> unique
        else
            varnames = [varname]
        end

        for varname in varnames
            if !haskey(plot_data, varname)
                plot_data[varname] = Dict{Symbol, Vector{Float64}}()
            end

            var_data = get!(plot_data[varname], structmod, Float64[])

            for nodesystem in nodesystems
                
                var_values = get_values(sol, nodesystem, varname)

                if var_values[1] isa Vector # multidimensional variable (e.g. vector of dimensions D)
                    for var_dimension in eachrow(reduce(hcat, var_values))
                        append!(var_data, var_dimension, NaN)
                    end
                else
                    append!(var_data, var_values, NaN)
                end

            end
        end
    end

    plots = []

    for varname in keys(plot_data)
        var_data = plot_data[varname] |> values |> collect
        group_sizes = length.(var_data)
        groups = [
            fill(structmod, group_size)
            for (structmod, group_size) in zip(keys(plot_data[varname]), group_sizes)
        ] |> x -> reduce(vcat, x)

        ys = reduce(vcat, var_data)
        xs = repeat(indep_values, length(ys) ÷ length(indep_values))
        push!(plots, plot(xs, ys, group = groups, title = "$varname"; kwargs...))
    end

    return length(plots) == 1 ? only(plots) : plots
end