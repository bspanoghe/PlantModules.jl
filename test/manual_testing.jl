using Revise, Infiltrator #!
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, MultiScaleTreeGraph
import PlantGraphs: Node
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots
using BenchmarkTools


begin # read in plant
	plant_graph = readXEG("./tutorials/temp/structures/beech10.xeg")
	mtg = convert_to_MTG(plant_graph)

	DataFrame(mtg, [:diameter, :length, :width])

	function combine_dimensions(l, d, w)
		if all(isnothing.([l, d, w]))
			return nothing
		elseif isnothing(w) # no width defined => shape is a cylinder
			return 1e2*[d/2, l] # 1e2 to go from m to cm
		else # otherwise we're dealing with a cuboid
			return 1e2*[l, w, 5e-4] # leaf thickness assumed to be 0.5 mm
		end
	end

	transform!(mtg, [:length, :diameter, :width] => combine_dimensions => :D)
	DataFrame(mtg, [:D]) # inspect the results

	traverse!(mtg, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")
	traverse!(mtg, node -> symbol!(node, "Stem"), symbol = "Internode")
	mtg = delete_nodes!(mtg, filter_fun = node -> isnothing(node.D))

	insert_parent!(mtg, NodeMTG("<", "Root", 0, 0))
	mtg = MultiScaleTreeGraph.get_root(mtg)

	if length(mtg) > 100 #!
		toomuch = [node for node in PlantModules.getnodes(mtg[1][1][1])]
		mtg = delete_nodes!(mtg, filter_fun = node -> node in toomuch)
	end
	DataFrame(mtg)
	plant_graph = mtg
end

if !isdefined(Main, :plant_graph) # manual plant
	mutable struct Root <: Node end
	mutable struct Stem <: Node end
	mutable struct Leaf <: Node
		D::Vector
	end

	plant_graph = Root() + Stem() + (Leaf([2]), Leaf([2.5]))	
end

# draw(plant_graph, resolution = (500, 400))
struct Soil <: Node end
struct Air <: Node end

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = PlantStructure(graphs, intergraph_connections)

default_changes = Dict(:Γ => 0.4, :T => 293.15)

module_defaults = Dict(
	:Root => Dict(:M => 300e-6),
	:Stem => Dict(:shape => Cylinder([6.0, 15.0], [0.8, 0.03]), :M => 400e-6),
	:Shoot => Dict(:shape => Cylinder([3.0, 8.0], [0.3, 0.01]), :M => 400e-6),
	:Leaf => Dict(:shape => Cuboid([5.0, 5.0, 20.0], [0.02, 0.02, 0.01]), :M => 450e-6, :K_s => 1e-5),
	:Soil => Dict(:W_max => 500.0, :T => 288.15),
	:Air => Dict(:W_r => 0.7, :K => 1e-4)
)

connecting_modules = [
	(:Soil, :Root) => (const_hydraulic_connection, Dict(:K => 500)),
	(:Root, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Leaf) => (const_hydraulic_connection, Dict()),
	(:Stem, :Leaf) => (const_hydraulic_connection, Dict(:K => 100)),
	(:Leaf, :Air) => (evaporation_connection, Dict()),
	(:Soil, :Air) => (hydraulic_connection, Dict())
]

func_connections = PlantFunctionality(; default_changes, module_defaults, connecting_modules)

module_coupling = Dict(
	:Root => [hydraulic_module, constant_carbon_module, K_module],
	:Stem => [hydraulic_module, constant_carbon_module, K_module],
	:Shoot => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
)

system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, [], (0.0, 5*24), sparse = true)
@time sol = solve(prob);

plotgraph(sol, graphs[1], varname = :W)
plotgraph(sol, graphs[1:2], varname = :Ψ)
plotgraph(sol, graphs[1], varname = :K, structmod = :Leaf)

plotgraph(sol, graphs[1], varname = [:D, :V], structmod = [:Stem, :Leaf]) .|> display
plotnode(sol, PlantModules.getnodes(graphs[1])[1], varname = [:P, :V]) .|> display

plotgraph(sol, graphs[1]) .|> display
plotgraph(sol, graphs[2]) .|> display
plotgraph(sol, graphs[3]) .|> display

# environmental_module(; name = :a, T = 280, W_max = 200.0, W_r = 0.5)
# environmental_module(; name = :b, T = 280, W_max = 99999999, W_r = 0.5)