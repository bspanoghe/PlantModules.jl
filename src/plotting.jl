# Plot MTK solutions #

"""
    plotgraph(sol::ODESolution, graph; structmod, varname, kwargs...)

Returns a plot for every functional variable for every node of the given graph for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable,
give the name of a structural module to limit considered nodes to those of this type, or both.
"""
function plotgraph(sol::ODESolution, graph; structmod = missing, varname = missing, kwargs...)
    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
    append!(indep_values, NaN) # NaN used to cause linebreaks in plot
    
    graphnodes = get_graphnodes(graph)
    node_structmods = PlantModules.getstructmod.(graphnodes)
    if !ismissing(structmod)
        node_structmods, graphnodes = filter_structmods(structmod, node_structmods, graphnodes)
    end

    nodesystems = getnodesystem.([sol], graphnodes)
    varname_dict = get_varname_dict(node_structmods, nodesystems, varname)
    varlist = [
        getproperty(nodesystems[nidx], _varname)
        for nidx in eachindex(nodesystems) for _varname in varname_dict[node_structmods[nidx]]
    ] # actual MTK variable names
    isempty(varlist) && error("Variable $varname not found in graph.")
    
    varlocs = getvarlocs(node_structmods, varname_dict, varlist)
    varvalues = sol[reduce(vcat, varlist)] |> x -> reduce(hcat, x) |> x -> [x fill(NaN, size(x, 1))]
    
    plots = AbstractPlot[]
    
    for _varname in unique(vcat(values(varname_dict)...))
        curr_varlocs = [varlocs[_structmod][_varname] for _structmod in keys(varlocs)] # vector per structmod with indexes of var values
        
        ys = varvalues[vcat(curr_varlocs...), :]' |> x -> vcat(x...)
        xs = repeat(indep_values, length(ys) ÷ length(indep_values))
        groups = [fill(_structmod, length(indep_values)*group_size) for (_structmod, group_size) in zip(keys(varlocs), length.(curr_varlocs))] |> x -> vcat(x...)
    
        push!(plots, plantplot(xs, ys, group = groups, title = "$_varname"; kwargs...))
    end

    return length(plots) == 1 ? only(plots) : plots
end

"""
    plotnode(sol::ODESolution, node; varname::Symbol)

Returns a plot for every functional variable of the given node for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable.
"""
function plotnode(sol::ODESolution, node; varname = missing, kwargs...)
    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
    append!(indep_values, NaN) # NaN used to cause linebreaks in plot

    structmod = PlantModules.getstructmod(node)
    nodesystem = getnodesystem(sol, node)

    varname_dict = get_varname_dict([structmod], [nodesystem], varname)
    varlist = [getproperty(nodesystem, _varname) for _varname in varname_dict[structmod]] # actual MTK variable names
    isempty(varlist) && error("Variable $varname not found in graph.")

    varlocs = getvarlocs([structmod], varname_dict, varlist)
    varvalues = sol[reduce(vcat, varlist)] |> x -> reduce(hcat, x) |> x -> [x fill(NaN, size(x, 1))]

    plots = Vector{AbstractPlot}(undef, length(varname_dict[structmod]))

    for (i, _varname) in enumerate(varname_dict[structmod])
        curr_varlocs = [varlocs[structmod][_varname]] # vector per structmod with indexes of var values
        
        ys = varvalues[vcat(curr_varlocs...), :]' |> x -> vcat(x...)
        xs = repeat(indep_values, length(ys) ÷ length(indep_values))
    
        plots[i] = plantplot(xs, ys, title = "$_varname"; kwargs...)
    end
  
    return length(plots) == 1 ? only(plots) : plots
end

# get all nodes of one or more graphs
get_graphnodes(graph) = PlantModules.getnodes(graph)
get_graphnodes(graph::Vector) = reduce(vcat, get_graphnodes.(graph))

# filter nodes of graph according to the structural module specified by the user
function filter_structmods(structmod::Symbol, node_structmods, graphnodes)
    chosen_structmods = node_structmods .== structmod
    if !any(chosen_structmods)
        error("Structural module \"$(structmod)\" not found in graph.")
    end

    return node_structmods[chosen_structmods], graphnodes[chosen_structmods]
end
function filter_structmods(structmod::Vector{Symbol}, node_structmods, graphnodes)
    chosen_structmods = [node_structmod in structmod for node_structmod in node_structmods]
    if !any(chosen_structmods)
        error("None of the structural modules \"$(structmod)\" were found in the graph.")
    end

    return node_structmods[chosen_structmods], graphnodes[chosen_structmods]
end

## Pry the ODE system corresponding with given node out of the ODE solution
function getnodesystem(sol::ODESolution, node)
    system = getfield(sol.prob.f, :sys)
    nodename = string(PlantModules.getstructmod(node)) * string(PlantModules.getid(node))
    nodesystem = [subsys for subsys in getproperty(system, :systems) if get_name(subsys) == Symbol(nodename)][1]
    return nodesystem
end

# returns what variables should be plotted per structmod
# e.g.: :Stem => [:W, :P, :M]
function get_varname_dict(node_structmods, nodesystems, varname)
    varname_dict = Dict{Symbol, Vector{Symbol}}() 
        
    for structmod_idx in unique(i -> node_structmods[i], eachindex(node_structmods))
        structmod = node_structmods[structmod_idx]
        varname_dict[structmod] = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[structmod_idx])] |> unique
        filter_varname_dict!(varname_dict, varname, structmod)
    end

    return varname_dict
end

# Get the symbol representation of a MTK unknown (variable)
get_MTKunknown_symbol(s::SymbolicUtils.Symbolic) = s.metadata[ModelingToolkit.VariableSource][2]

# filter varname_dict so only variable names specified by user remain, e.g. `Stem => [:W, :V]` => `Stem => [:V]`
filter_varname_dict!(varname_dict, varname::Missing, structmod) = nothing
function filter_varname_dict!(varname_dict, varname::Symbol, structmod)
    varname_dict[structmod] = [vn for vn in varname_dict[structmod] if vn == varname]
end
function filter_varname_dict!(varname_dict, varnames::Vector{Symbol}, structmod)
    for varname in varnames # user specified e.g. varname = :W
        varname_dict[structmod] = [vn for vn in varname_dict[structmod] if vn in varnames]
    end
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

# use RecipesBase.jl to avoid Plots.jl dependency
@userplot PlantPlot
@recipe function f(pp::PlantPlot)
    xs = pp.args[1]
    ys = pp.args[2]
    @series begin
        (xs, ys)
    end
end