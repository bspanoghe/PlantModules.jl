# # Plot plant structure

"""
    plotstructure(structure::PlantStructure)

Visualise the structure of a plant system.

This function is a thin wrapper around [`GraphRecipes.graphplot`](@ref) with some keyword arguments specified for plant structures.
"""
function plotstructure(structure::PlantStructure; kwargs...)
    names = getstructmod.(getnodes(structure))
    colordict = [name => idx for (idx, name) in enumerate(unique(names))] |> Dict
    markercolor = [colordict[name] for name in names]
    curves = false
    
    graphplot(structure; names, markercolor, curves, kwargs...)
end

# # Plot MTK solutions

"""
    plotgraph(sol::ODESolution, graph; structmod, varname, kwargs...)

Return a plot for every functional variable for every node of the given graph for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable,
give the name of a structural module to limit considered nodes to those of this type, or both.
"""
function plotgraph(sol::ODESolution, graph; structmod = missing, varname = missing, kwargs...)
    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
    append!(indep_values, NaN) # NaN used to cause linebreaks in plot
    
    graphnodes = getnodes(graph)
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
        for _structmod in keys(varlocs)
            @assert _varname in keys(varlocs[_structmod]) "$(_structmod) does not have the variable $(_varname) defined."
        end
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
    nodename = string(PlantModules.getstructmod(node)) * string(PlantModules.getid(node))

    system = sol.prob.f.sys
    parentsystem = get_parent(system) # system before simplification
    nodesystem = [subsys for subsys in get_systems(parentsystem) if get_name(subsys) == Symbol(nodename)][1]
    
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

#TODO: make own graph drawing recipe that works better for visualising tree graphs with cycles added at the ends

# @userplot StructurePlot

# function get_adj_matrix(ps::PlantStructure)
#     vs = PlantModules.vertices(ps)
#     n = length(vs)
#     adj_matrix = zeros(Bool, n, n)
#     for vertex in vs
#         neighbors = PlantModules.neighbors(ps, vertex)
#         adj_matrix[vertex, neighbors] .= true
#         adj_matrix[neighbors,  vertex] .= true
#     end
#     return adj_matrix
# end

# function get_weight_matrix(ps::PlantStructure)
#     vertices = PlantModules.vertices(ps)
#     num_neighbors = [length(PlantModules.neighbors(ps, vertex)) for vertex in vertices]

#     weight_matrix = [
#         1 / max(num_neighbors[i], num_neighbors[j])
#         for i in vertices, j in vertices
#     ]

#     return weight_matrix
# end

# function get_edge_positions(positions, ps::PlantStructure)
#     position_xs = first.(positions)
#     position_ys = last.(positions)

#     edge_xs = []
#     edge_ys = []

#     for vertex in PlantModules.vertices(ps)
#         for neighbor in PlantModules.neighbors(ps, vertex)
#             push!(edge_xs, [position_xs[vertex], position_xs[neighbor], missing])
#             push!(edge_ys, [position_ys[vertex], position_ys[neighbor], missing])
#         end
#     end

#     return (edge_xs, edge_ys)
# end

# @recipe function f(sp::StructurePlot)
#     plantsystem = sp.args[1]
#     adj_matrix = get_adj_matrix(plantsystem)
#     weight_matrix = get_weight_matrix(plantsystem)    
#     positions = GraphRecipes.NetworkLayout.stress(adj_matrix, weights = weight_matrix)
#     edge_positions = get_edge_positions(positions, plantsystem)

#     names = PlantModules.getstructmod.(PlantModules.getnodes(plantsystem))
#     colordict = [name => idx for (idx, name) in enumerate(unique(names))] |> Dict
#     markercolor = [colordict[name] for name in names]

#     @series begin
#         seriestype := :path
#         linecolor := :black
#         label := false
#         edge_positions
#     end

#     @series begin
#         seriestype := :scatter
#         label := false
#         markercolor := markercolor
#         markersize := 8
#         markershape := :hexagon
#         (first.(positions), last.(positions))
#     end
# end