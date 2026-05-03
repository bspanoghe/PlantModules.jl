# # Plot plant structure

"""
    plotstructure(plantstructure::PlantStructure; kwargs...)

Visualise the structure of a plant system.
"""
function plotstructure(plantstructure::PlantStructure; kwargs...)
    return structureplot(plantstructure; kwargs...)
end

plotstructure(graph; kwargs...) = plotstructure(PlantStructure(graph); kwargs...)


@userplot StructurePlot

function get_adj_matrix(ps::PlantStructure)
    vertices = getnodes(ps)

    adj_matrix = [
        vertices[i] in getneighbors(vertices[j], ps)
        for i in eachindex(vertices), j in eachindex(vertices)
    ]

    return adj_matrix
end

function get_edge_positions(positions, ps::PlantStructure)
    position_xs = first.(positions)
    position_ys = last.(positions)

    edge_xs = []
    edge_ys = []

    for vertex in PlantModules.vertices(ps)
        for neighbor in PlantModules.neighbors(ps, vertex)
            push!(edge_xs, [position_xs[vertex], position_xs[neighbor], missing])
            push!(edge_ys, [position_ys[vertex], position_ys[neighbor], missing])
        end
    end

    return (edge_xs, edge_ys)
end

@recipe function f(sp::StructurePlot)
    plantsystem = sp.args[1]

    # calculate positions
    adj_matrix = get_adj_matrix(plantsystem)
    positions = NetworkLayout.stress(adj_matrix)
    edge_positions = get_edge_positions(positions, plantsystem)
    xs = first.(positions)
    ys = last.(positions)

    # group nodes per structural module
    names = getstructmod.(getnodes(plantsystem))
    colordict = [name => idx for (idx, name) in enumerate(unique(names))] |> Dict
    markercolor = [colordict[name] for name in names]

    # set global plot attributes
    xmin, xmax = extrema(xs)
    Δx = xmax - xmin
    ymin, ymax = extrema(ys)
    Δy = ymax - ymin

    xlims --> (xmin - 0.1*Δx, xmax + 0.1*Δx)
    ylims --> (ymin - 0.1*Δy, ymax + 0.1*Δy)

    # plot edges
    @series begin 
        seriestype := :path
        linecolor := :black
        label := false
        edge_positions
    end

    # plot nodes
    @series begin 
        seriestype := :scatter
        label := false
        markercolor --> markercolor
        markersize := 10
        markershape := :hexagon
        (xs, ys)
    end

    # create legend
    @series begin 
        seriestype := :scatter
        label := sort(unique(names)) .|> string |> permutedims
        markercolor := [colordict[name] for name in sort(unique(names))] |> permutedims
        markersize := 6
        markershape := :hexagon
        (fill(NaN, (1, length(unique(names)))))
    end
end


# # Plot MTK solutions

"""
    plotgraph(sol::ODESolution, plantstructure::PlantStructure, nodes::Vector = getnodes(plantstructure); varname = missing, structmod = missing, kwargs...)

Return a plot for every functional variable for a collection of nodes of a plantstructure for the given solution `sol`.

Optionally, the user can give the name of a functional variable to only return a plot of this variable,
give the name of a structural module to limit considered nodes to those of this type, or both.
Alternatively, the data for the plot can be acquired directly with the function [`getplotdata`](@ref).
"""
function plotgraph(sol::ODESolution, plantstructure::PlantStructure, nodes::Vector = getnodes(plantstructure); varname = missing, structmod = missing, kwargs...)
    
    if varname isa Vector
        return [plotgraph(sol, plantstructure, nodes; varname = _varname, structmod, kwargs...) for _varname in varname]
    else
        xs, ys, groups = getplotdata(sol, plantstructure, nodes; varname, structmod)
        return plantplot(xs, ys, group = groups, title = "$varname"; kwargs...)
    end

end

plotgraph(sol::ODESolution, plantstructure::PlantStructure, node; varname = missing, structmod = missing, kwargs...) = (
    plotgraph(sol, plantstructure, [node]; varname, structmod, kwargs...)
) 
plotgraph(sol::ODESolution, graph, nodes::Vector = getnodes(PlantStructure(graph)); varname = missing, structmod = missing, kwargs...) = (
    plotgraph(sol, PlantStructure(graph), nodes; varname, structmod, kwargs...)
)
plotgraph(sol::ODESolution, graph, node; varname = missing, structmod = missing, kwargs...) = (
    plotgraph(sol, PlantStructure(graph), [node]; varname, structmod, kwargs...)
)

"""
    plotgraph(sols::Vector{<:ODESolution}, plantstructures::Vector{<:PlantStructure}, nodes_vec::Vector{<:Vector} = [getnodes(plantstructure) for plantstructure in plantstructures];
        varname = missing, structmod = missing, kwargs...)

Plot solutions for a collection of plantstructures and solutions, intended for sequential runs of a growing plant structure.

Note that if you want to plot specific nodes, you must pass a vector where each element is this set of nodes for every plantstructure in `plantstructures`.
"""
function plotgraph(sols::Vector{<:ODESolution}, plantstructures::Vector{<:PlantStructure}, nodes_vec::Vector{<:Vector} = [getnodes(plantstructure) for plantstructure in plantstructures];
        varname = missing, structmod = missing, kwargs...)
    
    if varname isa Vector
        return [plotgraph(sols, plantstructures, nodes_vec; varname = _varname, structmod, kwargs...) for _varname in varname] 
    else
        data = [getplotdata(sol, plantstructure, nodes; varname, structmod) for (sol, plantstructure, nodes) in zip(sols, plantstructures, nodes_vec)]
        xs_vec, ys_vec, groups_vec = [getindex.(data, i) for i in 1:3]
        for i in eachindex(xs_vec)[2:end]
            xs_vec[i] = xs_vec[i] .+ sum(sols[j].prob.tspan[2] for j in 1:i-1) # add ending time of previous solutions to times of current solution
        end
        xs, ys, groups = [reduce(vcat, i) for i in [xs_vec, ys_vec, groups_vec]]

        return plantplot(xs, ys, group = groups, title = "$varname"; kwargs...)
    end

end

plotgraph(sols::Vector{<:ODESolution}, plantstructures::Vector{<:PlantStructure}, nodes::Vector; varname = missing, structmod = missing, kwargs...) = (
    plotgraph(sols, plantstructures, [[node] for node in nodes]; varname, structmod, kwargs...)
)

"""
    getplotdata(sol::ODESolution, plantstructure::PlantStructure; varname, structmod, nodes)

Get the x-values, y-values and groups required to plot the given solution. See [`graphplot`](@ref) for more information. 
"""
function getplotdata(sol::ODESolution, plantstructure::PlantStructure, nodes; varname, structmod)
    indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
    append!(indep_values, NaN) # NaN used to cause linebreaks in plot

    varlist, node_structmods, varname_dict = _getvariables(sol, plantstructure, varname, structmod, nodes)

    varlocs = getvarlocs(node_structmods, varname_dict, varlist) # e.g.: varlocs[:Stem][:W] => [10, 15, 16]
    varvalues = sol[reduce(vcat, varlist)] |> x -> reduce(hcat, x) |> x -> [x fill(NaN, size(x, 1))]

    for _structmod in keys(varlocs)
        @assert varname in keys(varlocs[_structmod]) "$(_structmod) does not have the variable $(varname) defined."
    end
    curr_varlocs = [varlocs[_structmod][varname] for _structmod in keys(varlocs)] # vector per structmod with indexes of var values

    ys = varvalues[vcat(curr_varlocs...), :]' |> x -> vcat(x...)
    xs = repeat(indep_values, length(ys) ÷ length(indep_values))
    groups = [fill(_structmod, length(indep_values) * group_size) for (_structmod, group_size) in zip(keys(varlocs), length.(curr_varlocs))] |> x -> vcat(x...)

    return xs, ys, groups
end

getplotdata(sol::ODESolution, plantstructure::PlantStructure; varname, structmod) = getplotdata(sol, plantstructure, getnodes(plantstructure); varname, structmod)

"""
    getvariables(sol::ODESolution, plantstructure::PlantStructure; varname = missing, structmod = missing)

Return the Numeric representation of one or more variables from a plantstructure, optionally filtered by structural module.
"""
function getvariables(sol::ODESolution, plantstructure::PlantStructure; varname = missing, structmod = missing)
    varlist, _, _ = _getvariables(sol, plantstructure, varname, structmod)
    return varlist
end

# internal version with more outputs than users need
function _getvariables(sol::ODESolution, plantstructure::PlantStructure, varname, structmod, graphnodes::Vector)
    node_structmods = PlantModules.getstructmod.(graphnodes)

    if !ismissing(structmod)
        node_structmods, graphnodes = filter_structmods(structmod, node_structmods, graphnodes)
    end
    nodesystems = getnodesystem.([sol], graphnodes, [plantstructure])
    varname_dict = get_varname_dict(node_structmods, nodesystems, varname)
    varlist = [
        getproperty(nodesystems[nidx], _varname)
            for nidx in eachindex(nodesystems) for _varname in varname_dict[node_structmods[nidx]]
    ]
    isempty(varlist) && error("Variable $varname not found in graph.")

    return varlist, node_structmods, varname_dict
end

_getvariables(sol::ODESolution, plantstructure::PlantStructure, varname, structmod) = _getvariables(sol, plantstructure, varname, structmod, getnodes(plantstructure))


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
function getnodesystem(sol::ODESolution, node, plantstructure::PlantStructure)
    nodename = getsysname(node, plantstructure)
    sys = sol.prob.f.sys
    nodesystem = getsubsystem(sys, nodename)

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
get_MTKunknown_symbol(s) = (
    operation(s) == getindex ? 
	Symbol(operation(arguments(s)[findfirst(x -> iscall(x), arguments(s))])) :
	Symbol(operation(s))
)

# filter varname_dict so only variable names specified by user remain, e.g. `Stem => [:W, :V]` => `Stem => [:V]`
filter_varname_dict!(varname_dict, varname::Missing, structmod) = nothing
function filter_varname_dict!(varname_dict, varname::Symbol, structmod)
    return varname_dict[structmod] = [vn for vn in varname_dict[structmod] if vn == varname]
end
function filter_varname_dict!(varname_dict, varnames::Vector{Symbol}, structmod)
    for varname in varnames # user specified e.g. varname = :W
        varname_dict[structmod] = [vn for vn in varname_dict[structmod] if vn in varnames]
    end
    return
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
            varidxs = collect((get(cumulvarlengths, nc - 1, 0) + 1):cumulvarlengths[nc])
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