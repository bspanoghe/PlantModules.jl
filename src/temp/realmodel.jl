using Pkg; Pkg.activate(".")
include("../PlantModules.jl"); using .PlantModules
using PlantGraphs, ModelingToolkit, DifferentialEquations, Unitful, Plots, MultiScaleTreeGraph
import GLMakie.draw

using MutableNamedTuples

struct PMNode
	structmod::Symbol
	attributes::MutableNamedTuple
end

abstract type MyPGNode <: PlantGraphs.Node end

struct Node{T} <: MyPGNode
	structmod::Symbol
    attributes::Dict
end

Node(type::Symbol, attributes) = Node{type, typeof(attributes)}(attributes)

a = Node(:leaf, ["HAH!", "HOOH"])
b = Node(:shoot, (horses = 3, poseidon = :LORD_POSEIDON))

r = a + b

n = PlantModules.nodes(r)[1]
PlantModules.id(n)

structmods = [
	[PlantModules.structmod(node), collect(keys(PlantModules.attributes(node)))...]
	for node in PlantModules.nodes(me_graph)
] |> x -> unique(vec -> vec[1], x)

define_node(:moo, :woo)
@check_inputs(structmods[1])

# Structural modules #

me_graph = readXEG("./src/temp/structures/beech0.xeg") #!
convert_to_PG(graph) |> draw

## Plant
plant_graph = readXEG("./src/temp/structures/beech0.xeg") #! change to beech
plotbranching(plant_graph)

mtg = convert_to_MSTG(plant_graph)

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
	println("There are $num_nodes nodes of type $me_structmod $( all(isnothing.(dimensions)) ? "(dimensions undefined)" : "")")
end

traverse!(mtg, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")

descendants(mtg, symbol = "Internode", self = true) |> DataFrame
descendants(mtg, symbol = "Shoot", self = true) |> DataFrame
descendants(mtg, symbol = "Leaf", self = true) |> DataFrame

## Environment

struct Soil <: Node end
struct Air <: Node end

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

## connections

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now
struct_connections = [graphs, intergraph_connections]


# Functional modules #

## New functional modules

function photosynthesis(x)
	println("BAGOOL!") #!
end

## Connect them to structure

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

connecting_modules = [
	() => PlantModules.hydraulic_connection,
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 8]),
	(:Root, :Internode) => (PlantModules.hydraulic_connection, [:K => 4]), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Internode, :Shoot) => (PlantModules.hydraulic_connection, [:K => 1e-2]),
	(:Internode, :Leaf) => (PlantModules.hydraulic_connection, [:K => 5e-2]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

func_connections = [connecting_modules, get_connection_eqs]

## Tweak parameters

C_stem = 300
C_shoot = 350
C_leaf = 400

module_defaults = (
	Internode = (shape = PlantModules.Cilinder(ϵ_D = [5.0, 0.3], ϕ_D = [0.7, 0.1]), M = C_stem),
	Shoot = (shape = PlantModules.Cilinder(ϵ_D = [5.0, 0.3], ϕ_D = [0.7, 0.1]), M = C_shoot),
	Leaf = (shape = PlantModules.Cuboid(ϵ_D = [0.5, 0.5, 0.01], ϕ_D = [0.5, 0.5, 0.01]), M = C_leaf),
	Soil = (W_max = 2000.0, T = 288.15),
	Air = ()
)

# Gettem #

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))	 #! Generate warning for missing variable defaults that aren't dummies
sol = solve(prob)

PlantModules.plotgraph(sol, graphs[1]) .|> display
PlantModules.plotgraph(sol, graphs[2]) .|> display
PlantModules.plotgraph(sol, graphs[3]) .|> display