using Pkg; Pkg.activate(".")
include("../PlantModules.jl")
using PlantGraphs, ModelingToolkit, DifferentialEquations, Unitful, Plots, MultiScaleTreeGraph
import GLMakie.draw

# Creating plant #

lines = readlines("./src/temp/structures/beech.xeg")
htmlmatch(element, line) = match(Regex("(?<=$(element)=\").+?(?=\")"), line).match

newgraph = Dict{Int, Dict}()

curr_id = 0

for line in lines
	line_element = split(line)[1]
	if line_element == "<node"
		id = htmlmatch("id", line) |> x -> parse(Int, x)
		structmod = htmlmatch("type", line)

		newgraph[id] = Dict(:nb => Int[], :at => Dict{Symbol, Any}(), :sm => structmod)
		curr_id = id
	elseif line_element == "<edge" 
		id = htmlmatch("src_id", line) |> x -> parse(Int, x)
		nb_id = htmlmatch("dest_id", line) |> x -> parse(Int, x)
		push!(newgraph[id][:nb], nb_id)
	elseif line_element == "<property" && occursin("value=", line)
		name = htmlmatch("name", line)
		val = htmlmatch("value", line)
		if !isnothing(match(r"^[0-9]+$", val))
			val = parse(Int, val)
		elseif !isnothing(match(r"^-?[0-9]*\.[0-9]+(E-[0-9]+)?$", val))
			val = parse(Float64, val)
		end
		newgraph[curr_id][:at][Symbol(name)] = val
	end
end

nodes(graph::Dict) = collect(graph)
neighbours(node::Pair, graph::Dict) = [graph[nbidx] for nbidx in node.second[:nb]]
attributes(node::Pair) = node.second[:at]
structmod(node::Pair) = node.second[:sm]
id(node::Pair) = node.first

newnode = nodes(newgraph)[11]
neighbours(newnode, newgraph)
attributes(newnode)
structmod(newnode)
id(newnode)

mutable struct Root <: Node 
	D::Vector
end

mutable struct PNode <: Node
end
mutable struct Internode <: Node
	D::Vector
end
mutable struct Meristem <: Node
	G
end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

## parameters
branching_prob = 0.1
G_th = 5
G0 = 10
rootsize0 = [0.05, 0.5]
internodesize0 = [0.05, 0.5]
leafsize0 = [0.5, 0.2, 0.01]

## rules

shoot_rule = Rule(Meristem, lhs = meristem -> rand() < data(meristem).G, rhs = meristem -> PNode() + (Meristem(branching_prob) + Leaf(leafsize0), Internode(internodesize0) + Meristem(1)))
# root_rule = Rule(Meristem, lhs = meristem -> rand() < data(meristem).G, rhs = meristem -> Root(rootsize0) + Meristem(branching_prob), Root(rootsize0) + Meristem(branching_prob), Root(rootsize0) + Meristem(branching_prob))

## grow 'em
axiom = PNode() + Internode(internodesize0) + Meristem(G0)

plant_graph = Graph(axiom = axiom, rules = (shoot_rule,))
[rewrite!(plant_graph) for _ in 1:10]
draw(plant_graph)


# E #

C_root = 250 # (mol/m^3) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 300 # Idem here
C_leaf = 330 # And here!
# loosely based on https://www.researchgate.net/figure/Total-soluble-sugar-content-in-the-leaf-a-and-root-tissues-b-of-Homjan-HJ_fig3_226164026

module_defaults = (
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.7, 0.1, 0.05]), D = [0.3, 0.03, 0.01], M = C_root),
	Pnode = (), #! it just aint right
    Internode = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.8, 0.03]), D = [0.015, 0.1], M = C_stem),
    Meristem = (),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 2000.0, T = 288.15),
	Air = ()
)

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]



soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = [graphs, intergraph_connections]

connecting_modules = [
	() => PlantModules.hydraulic_connection,
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 8]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 4]), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-2]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-2]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

func_connections = [connecting_modules, get_connection_eqs]

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))	 #! Generate warning for missing variable defaults that aren't dummies
sol = solve(prob)

PlantModules.plotgraph(sol, graphs[1]) .|> display
PlantModules.plotgraph(sol, graphs[2]) .|> display
PlantModules.plotgraph(sol, graphs[3]) .|> display