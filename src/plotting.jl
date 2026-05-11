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
    plotgraph(sol::ODESolution, plantstructure::PlantStructure, nodes::Vector = getnodes(plantstructure);
    varname::Union{Symbol, Missing} = missing, structmod::Union{Symbol, Vector{Symbol}, Missing} = missing, kwargs...
    )

Return a plot for a functional variable for a collection of nodes of a plantstructure for the given solution `sol`.

Optionally, the user can give the name of a structural module to limit considered nodes to those of this type.
Alternatively, the data for the plot can be acquired directly with the function [`getplotdata`](@ref).
"""
function plotgraph(sol::ODESolution, plantstructure::PlantStructure, nodes::Vector = getnodes(plantstructure);
    varname::Symbol, structmod::Union{Symbol, Vector{Symbol}, Missing} = missing, kwargs...)
    
    xs, ys, groups = getplotdata(sol, plantstructure, varname, filter_nodes(nodes, structmod))
    return plantplot(xs, ys; groups, title = "$varname", kwargs...)
end

"""
    getplotdata(sol::ODESolution, plantstructure::PlantStructure, nodes::Vector = getnodes(plantstructure); 
    varname::Symbol, structmod::Union{Symbol, Vector{Symbol}, Missing} = missing)

Get the x-values, y-values and groups required to plot the given solution. See [`graphplot`](@ref) for more information. 
"""
function getplotdata(sol::ODESolution, plantstructure::PlantStructure, nodes::Vector = getnodes(plantstructure); 
    varname::Symbol, structmod::Union{Symbol, Vector{Symbol}, Missing} = missing)

    nodes = filter_nodes(nodes, structmod)
    nodesystems = [getnodesystem(sol, node, plantstructure) for node in nodes]
    node_vars = [getproperty(nodesystem, varname) for nodesystem in nodesystems]
    node_values = sol[node_vars]

    xs = [copy(sol[get_iv(sol.prob.f.sys)]); NaN] |> # values of indepedent variable (NaN used to cause linebreaks)
        x -> repeat(x, length(nodes))
    ys = node_values' |> 
        x -> reduce(vcat, x) |>
        x -> [x; fill(NaN, 1, size(x, 2))] |> # add NaNs to ys for linebreaks as well
        x -> reduce(vcat, x)
    groups = [fill(getstructmod(node), length(node_values) + 1) for node in nodes] |> 
        x -> reduce(vcat, x)

    return xs, ys, groups
end

# filter nodes based on structural module
filter_nodes(nodes::Vector, ::Missing) = nodes
filter_nodes(nodes::Vector, structmod::Symbol) = [node for node in nodes if getstructmod(node) == structmod]
filter_nodes(nodes::Vector, structmods::Vector{Symbol}) = [node for node in nodes if getstructmod(node) in structmods]

# Pry the ODE system corresponding with given node out of the ODE solution
function getnodesystem(sol::ODESolution, node, plantstructure::PlantStructure)
    nodename = getsysname(node, plantstructure)
    sys = sol.prob.f.sys
    nodesystem = getsubsystem(sys, nodename)

    return nodesystem
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