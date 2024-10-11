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

get_graphnodes(graph) = PlantModules.nodes(graph)
get_graphnodes(graph::Vector) = reduce(vcat, get_graphnodes.(graph))

# what variables should be plotted per structmod
# e.g.: :Stem => [:W, :P, :M]
function getvarnames(structmods, nodesystems, varname)
    varnames = Dict{Symbol, Vector{Symbol}}() 
        
    for structmod_idx in unique(i -> structmods[i], eachindex(structmods))
        structmod = structmods[structmod_idx]
        varnames[structmod] = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[structmod_idx])] |> unique
        if varname != Symbol("") # user specified e.g. varname = :W
            varnames[structmod] = (varname in varnames[structmod]) ? [varname] : Symbol[]
        end
    end

    return varnames
end

# get rows of varlist that correspond with given structmod and varname
# e.g.: varlocs[:Stem][:W] => [10, 15, 16]
function getvarlocs(structmods, varnames, varlist)
    cumulvarlengths = length.(varlist) |> cumsum

    varlocs = Dict{Symbol, Dict{Symbol, Vector{Int64}}}()
    nc = 0 # nodecounter
    for _structmod in structmods
        if !haskey(varlocs, _structmod)
            varlocs[_structmod] = Dict{Symbol, Vector{Int64}}()
        end
    
        for _varname in varnames[_structmod]
            structlocs = get!(varlocs[_structmod], _varname, Int64[])
    
            nc += 1
            varidxs = collect(get(cumulvarlengths, nc-1, 0)+1:cumulvarlengths[nc])
            append!(structlocs, varidxs)
        end
    end

    return varlocs
end

"""
    plotnode(sol::ODESolution, node; varname::Symbol)

Returns a plot for every functional variable of the given node for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable.
"""
function plotnode(sol::ODESolution, node; varname::Symbol = Symbol(""), kwargs...)
    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
    append!(indep_values, NaN) # NaN used to cause linebreaks in plot

    structmod = PlantModules.structmod(node)
    nodesystem = getnodesystem(sol, node)

    varnames = getvarnames([structmod], [nodesystem], varname)
    varlist = [getproperty(nodesystem, _varname) for _varname in varnames[structmod]]
    varlocs = getvarlocs([structmod], varnames, varlist)
    varvalues = sol[reduce(vcat, varlist)] |> x -> reduce(hcat, x) |> x -> [x fill(NaN, size(x, 1))]

    plots = Vector{Plots.Plot}(undef, length(varnames[structmod]))

    for (i, _varname) in enumerate(varnames[structmod])
        curr_varlocs = [varlocs[structmod][_varname]] # vector per structmod with indexes of var values
        
        ys = varvalues[vcat(curr_varlocs...), :]' |> x -> vcat(x...)
        xs = repeat(indep_values, length(ys) ÷ length(indep_values))
    
        plots[i] = plot(xs, ys, title = "$_varname"; kwargs...)
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

    graphnodes = get_graphnodes(graph)
    
    structmods = PlantModules.structmod.(graphnodes)
    if structmod != Symbol("") # user specified e.g. `structmod = :Leaf` => filter out other nodes
        chosen_structmods = structmods .== structmod
        structmods = structmods[chosen_structmods]
        graphnodes = graphnodes[chosen_structmods]
    end
    
    nodesystems = getnodesystem.([sol], graphnodes)
    varnames = getvarnames(structmods, nodesystems, varname)
    
    varlist = [
        getproperty(nodesystems[nidx], _varname)
        for nidx in eachindex(nodesystems) for _varname in varnames[structmods[nidx]]
    ]
    
    varlocs = getvarlocs(structmods, varnames, varlist)
    varvalues = sol[reduce(vcat, varlist)] |> x -> reduce(hcat, x) |> x -> [x fill(NaN, size(x, 1))]
    
    plots = Plots.Plot[]
    
    for _varname in unique(vcat(values(varnames)...))
        curr_varlocs = [varlocs[_structmod][_varname] for _structmod in keys(varlocs)] # vector per structmod with indexes of var values
        
        ys = varvalues[vcat(curr_varlocs...), :]' |> x -> vcat(x...)
        xs = repeat(indep_values, length(ys) ÷ length(indep_values))
        groups = [fill(_structmod, length(indep_values)*group_size) for (_structmod, group_size) in zip(keys(varlocs), length.(curr_varlocs))] |> x -> vcat(x...)
    
        push!(plots, plot(xs, ys, group = groups, title = "$_varname"; kwargs...))
    end

    return length(plots) == 1 ? only(plots) : plots
end