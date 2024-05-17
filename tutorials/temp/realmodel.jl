using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs, ModelingToolkit, DifferentialEquations, Unitful, Plots, MultiScaleTreeGraph
import GLMakie.draw

# Structural modules #

import MultiScaleTreeGraph: delete_nodes!, delete_nodes_!

# function delete_node!(node::Node{N,A}; child_link_fun=new_child_link) where {N<:AbstractNodeMTG,A}
#     if isroot(node)
#         if length(children(node)) == 1
#             # If it has only one child, make it the new root:
#             chnode = children(node)[1]
#             # Add to the new root the mandatory root attributes:
#             root_attrs = Dict(
#                 :symbols => node[:symbols],
#                 :scales => node[:scales],
#                 :description => node[:description]
#             )

#             append!(chnode, root_attrs)

#             link!(chnode, child_link_fun(chnode))
#             reparent!(chnode, nothing)

#             node_return = chnode
#         else
#             error("Can't delete the root node if it has several children")
#         end
#     else
#         parent_node = parent(node)

#         if !isleaf(node)
#             # We re-parent the children to the parent of the node.
#             for chnode in children(node)
#                 # Updating the link of the children:
#                 link!(chnode, child_link_fun(chnode))
#                 addchild!(parent_node, chnode; force=true)
#             end
#         end

#         # Delete the node as child of his parent:
#         deleteat!(children(parent_node), findfirst(x -> node_id(x) == node_id(node), children(parent_node)))
#         node_return = parent_node
#     end

#     node = nothing

#     return node_return
# end

function delete_nodes!(
    node;
    scale=nothing,
    symbol=nothing,
    link=nothing,
    all::Bool=true, # like continue in the R package, but actually the opposite
    filter_fun=nothing,
    child_link_fun=new_child_link
	)

    # Check the filters once, and then compute the descendants recursively using `descendants_`
    check_filters(node, scale=scale, symbol=symbol, link=link)
    filtered = is_filtered(node, scale, symbol, link, filter_fun)

    while filtered
        node = delete_node!(node)
        filtered = is_filtered(node, scale, symbol, link, filter_fun)

        # Don't go further if all == false
        !all && return
    end

    delete_nodes!_(node, scale, symbol, link, all, filter_fun, child_link_fun)

    return node
end

function delete_nodes!_(node, scale, symbol, link, all, filter_fun, child_link_fun)
    if !isleaf(node)
        # First we apply the algorithm recursively on the children:
        chnodes = children(node)
        nchildren = length(chnodes)
        #? Note: we don't use `for chnode in chnodes` because it may delete dynamically during traversal, so we forget to traverse some nodes
        for chnode in chnodes[1:nchildren]
            delete_nodes!_(chnode, scale, symbol, link, all, filter_fun, child_link_fun)
        end
    end

    # Then we work on the node itself. This ensures that its children will not be deleted
    # afterwards (the deletion is acropetal, i.e. from leaves to root)

    # Is there any filter happening for the current node? (true is deleted):
    filtered = is_filtered(node, scale, symbol, link, filter_fun)

    if filtered
        delete_node!(node, child_link_fun=child_link_fun)
    end
end


## Plant
plant_graph = readXEG("tutorials/temp/structures/beech3.xeg") #! change to beech
# convert_to_PG(plant_graph) |> draw

mtg = convert_to_MTG(plant_graph)

### Setting attributes right
DataFrame(mtg, [:diameter, :length, :width])

function combine_dimensions(l, d, w)
	if all(isnothing.([l, d, w]))
		return nothing
	elseif isnothing(w)
		return [l, d]
	else
		return [l, w, 1.5e-4]
	end
end

transform!(mtg, [:length, :diameter, :width] => combine_dimensions => :D)
DataFrame(mtg, [:D])

### Inspecting what kind of structural modules are in here
me_structmods = [PlantModules.structmod(node) for node in PlantModules.nodes(mtg)] |> unique

for me_structmod in me_structmods
	dimensions = [node_attributes(node)[:D] for node in PlantModules.nodes(mtg) if PlantModules.structmod(node) == me_structmod]
	num_nodes = length(dimensions)
	println("There are $num_nodes nodes of type $me_structmod $( all(ismissing.(dimensions)) ? "(dimensions undefined)" : "")")
end

traverse!(mtg, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")
mtg = delete_nodes!(mtg, filter_fun = node -> isnothing(node.D))

mtg |> DataFrame

descendants(mtg, symbol = "Internode", self = true) |> DataFrame
descendants(mtg, symbol = "Shoot", self = true) |> DataFrame
descendants(mtg, symbol = "Leaf", self = true) |> DataFrame

## Environment

struct Soil <: PlantGraphs.Node end
struct Air <: PlantGraphs.Node end

soil_graph = Soil()
air_graph = Air()

graphs = [mtg, soil_graph, air_graph]

## connections

intergraph_connections = [[1, 2] => (mtg, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)]
struct_connections = [graphs, intergraph_connections]

# Functional modules #

## New functional modules

function photosynthesis(x)
	println("BAGOOL!") #!
end

## Connect them to structure

module_coupling = [
	PlantModules.hydraulic_module => [:Internode, :Shoot, :Leaf],
	PlantModules.constant_carbon_module => [:Internode, :Shoot, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

connecting_modules = [
	() => PlantModules.hydraulic_connection,
	(:Soil, :Internode) => (PlantModules.hydraulic_connection, [:K => 100]),
    (:Internode, :Internode) => (PlantModules.hydraulic_connection, [:K => 3]),
	(:Internode, :Shoot) => (PlantModules.hydraulic_connection, [:K => 2]),
	(:Internode, :Leaf) => (PlantModules.hydraulic_connection, [:K => 2]),
    (:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 0])
]

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

func_connections = [connecting_modules, get_connection_eqs]

## Tweak parameters

C_stem = 300
C_shoot = 350
C_leaf = 400

module_defaults = (
	Internode = (shape = PlantModules.Cilinder(ϵ_D = [5.0, 0.3], ϕ_D = [0.1, 0.01]), M = C_stem),
	Shoot = (shape = PlantModules.Cilinder(ϵ_D = [5.0, 0.3], ϕ_D = [0.1, 0.01]), M = C_shoot),
	Leaf = (shape = PlantModules.Cuboid(ϵ_D = [0.5, 0.5, 0.01], ϕ_D = [0.1, 0.1, 0.01]), M = C_leaf),
	Soil = (W_max = 100000.0, T = 288.15),
	Air = ()
)

# Gettem #

system = PlantModules.generate_system(PlantModules.default_params, PlantModules.default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system)
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
sol = solve(prob)

PlantModules.plotgraph(sol, graphs[1], func_varname = :W)
PlantModules.plotgraph(sol, graphs[2], func_varname = :W)
PlantModules.plotgraph(sol, graphs[1], func_varname = :Ψ)

# plotspeed

xs = [rand(100) for _ in 1:1000]
ys = [rand(100) for _ in 1:1000]

function plot1(xs, ys)
    plt = plot()
    for (x, y) in zip(xs, ys)
        plot!(plt, x, y)
    end
    plt
end

function plot2(xs, ys)
    plot(reduce(vcat, xs), reduce(vcat, ys))
end

@btime plot1($xs, $ys) ; # 431.335 ms (2231129 allocations: 121.87 MiB)
@btime plot2($xs, $ys) ; # 2.844 ms (2530 allocations: 3.20 MiB)

xs = [i%2 == 0 ? NaN : 1:10 for i in 1:6]
ys = [i%2 == 0 ? NaN : rand(10) for i in 1:6]
plot2(xs, ys)