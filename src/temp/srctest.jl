using Pkg; Pkg.activate(".")
include("../PlantModules.jl")
using PlantGraphs, ModelingToolkit, DifferentialEquations, Unitful, Plots #! MTK imports etc. should not be necessary when package is done
# import GLMakie.draw
# using BenchmarkTools

#! for testing
function showme(x)
	[println(xi) for xi in sort(collect(string.(x)))]
	return nothing
end

mutable struct Root <: Node end
mutable struct Stem <: Node end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

C_root = 250 # (mol/m^3) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 300 # Idem here
C_leaf = 330 # And here!
# loosely based on https://www.researchgate.net/figure/Total-soluble-sugar-content-in-the-leaf-a-and-root-tissues-b-of-Homjan-HJ_fig3_226164026

PlantModules.default_params
PlantModules.default_u0s

model_defaults = (Γ = 0.4, P = 0.2, T = 293.15) #! currently assumed variables with same name are the same for all functional modules
#! replace with function that generates new default_params

module_defaults = (
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.7, 0.1, 0.05]), D = [0.3, 0.03, 0.01], M = C_root), #! changed C to M_const #! include check whether module_defaults variabkle names correspond with default_params / default_u0s ?
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.8, 0.03]), D = [0.015, 0.1], M = C_stem),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 500.0, T = 288.15),
	Air = ()
)

module_coupling = [ #! switch around?
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

if !isdefined(Main, :plant_graph) # for debugging mostly
	plant_graph = Root() + Stem() + (Leaf([0.02]), Leaf([0.025]))
end

# draw(plant_graph, resolution = (500, 400))

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

system = PlantModules.generate_system(model_defaults, module_defaults, module_coupling,
	struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))	 #! Generate warning for missing variable defaults that aren't dummies
sol = solve(prob)

plot(sol, idxs = [Symbol("Soil1₊Ψ"), Symbol("Root1₊Ψ"), Symbol("Stem2₊Ψ"), Symbol("Leaf3₊Ψ"), Symbol("Leaf4₊Ψ"), Symbol("Air1₊Ψ")])
plot(sol, idxs = [Symbol("Soil1₊Ψ"), Symbol("Root1₊Ψ"), Symbol("Stem2₊Ψ"), Symbol("Leaf3₊Ψ"), Symbol("Leaf4₊Ψ")])
plot(sol, idxs = [Symbol("Root1₊Π"), Symbol("Stem2₊Π"), Symbol("Leaf3₊Π"), Symbol("Leaf4₊Π")])
plot(sol, idxs = [Symbol("Root1₊P"), Symbol("Stem2₊P"), Symbol("Leaf3₊P"), Symbol("Leaf4₊P")])

plot(sol, idxs = [Symbol("Root1₊M")])
plot(sol, idxs = [Symbol("Leaf3₊M")])
plot(sol, idxs = [Symbol("Leaf4₊M")])

plot(sol, idxs = [Symbol("Soil1₊W")])
plot(sol, idxs = [Symbol("Root1₊W")])
plot(sol, idxs = [Symbol("Stem2₊W")])
plot(sol, idxs = [Symbol("Leaf3₊W")])
plot(sol, idxs = [Symbol("Leaf4₊W")])
plot(sol, idxs = [Symbol("Air1₊W")])

plot(sol, idxs = [Symbol("Root1_Soil1₊F"), Symbol("Air1_Leaf3₊F")])

# tests #


# Plotting tests #

function lotka_volterra(; name, α, β, δ, γ, N, P)
	@parameters α = α β = β δ = δ γ = γ
	@variables t N(t) = N P(t) = P ΣF_N(t) ΣF_P(t)
	d = Differential(t);
	eqs = [
		d(N) ~ α*N - β*N*P + ΣF_N,
		d(P) ~ δ*N*P - γ*P + ΣF_P
	]
	return ODESystem(eqs, t; name)
end

function fountain_of_rabbits(; name, η, N, P)
	@variables t N(t) = N P(t) = P ΣF_N(t) ΣF_P(t)
	d = Differential(t);
	eqs = [
		d(N) ~ η + ΣF_N,
		d(P) ~ ΣF_P
	]
	return ODESystem(eqs, t; name)
end

function wandering_animals(; name, κ)
    @parameters (
        κ = κ,
    )
    @variables (
		t,
        F_N(t),
        N_1(t),
        N_2(t),
		F_P(t),
		P_1(t),
		P_2(t)
    )

    eqs = [
        F_N ~ κ * (N_2 - N_1),
		F_P ~ κ * (P_2 - P_1)
    ]
    return ODESystem(eqs, t; name)
end

wandering_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
	[connection_MTK.N_1 ~ node_MTK.N for connection_MTK in connection_MTKs]...,
	[connection_MTK.P_1 ~ node_MTK.P for connection_MTK in connection_MTKs]...,
	[connection_MTK.N_2 ~ nb_node_MTK.N for (connection_MTK, nb_node_MTK) in zip(connection_MTKs, nb_node_MTKs)]...,
	[connection_MTK.P_2 ~ nb_node_MTK.P for (connection_MTK, nb_node_MTK) in zip(connection_MTKs, nb_node_MTKs)]...,
	node_MTK.ΣF_N ~ sum([connection_MTK.F_N for connection_MTK in connection_MTKs]),
	node_MTK.ΣF_P ~ sum([connection_MTK.F_P for connection_MTK in connection_MTKs])
]

Base.@kwdef mutable struct Forest <: Node
	N::Int = 20
	P::Int = 10
	δ::Float64 = 1.0
end

Base.@kwdef mutable struct Grassland <: Node
	N::Int = 100
	P::Int = 5
end

Base.@kwdef mutable struct Cave <: Node end

default_params = (
	lotka_volterra = (α = 1.5, β = 2.0, δ = 1.5, γ = 0.8),
	fountain_of_rabbits = (η = 10,),
	wandering_animals = (κ = 0.1,)
)
default_u0s = (
	lotka_volterra = (N = 30, P = 10),
	fountain_of_rabbits = (N = 1, P = 0),
	wandering_animals = ()
)
model_defaults = (β = 1.9, P = 15)
module_defaults = (
	Grassland = (β = 0.1,),
	Forest = (δ = 2.3,),
	Cave = (β = 0,)
)
module_coupling = [
	lotka_volterra => [:Grassland, :Forest],
	fountain_of_rabbits => [:Cave]
]
graph = Grassland() + Grassland(N = 200, P = 3) + (Forest(δ = 1.8), Forest(N = 5, P = 0, δ = 2.3))
graphs = [graph, Cave()]
struct_connections = [graphs, [[1, 2] => (:Forest, :Cave)]] #! make vector of graphs as input unnecessary when theres only 1 graph
func_connections = [[() => wandering_animals], wandering_eqs]

graphs, intergraph_connections = struct_connections
connecting_modules, get_connection_eqs = func_connections


# graphfuncs #

node1, node2, node3, node4 = collect(values(graph.nodes))

## nodes
issetequal(PlantModules.nodes(graph), [node1, node2, node3, node4])

## neighbours
issetequal(PlantModules.neighbours(node2, graph), [node1, node3, node4])

## attributes
issetequal(PlantModules.attributes(node4), [:P => 0, :N => 5, :δ => 2.3])

## structmod
PlantModules.structmod(node1) == :Grassland

## id
allunique(PlantModules.id.([node1, node2, node3, node4]))


# generate_system #

## getMTKsystem
node1, node2, node3, node4 = collect(values(graph.nodes))

sys1 = PlantModules.getMTKsystem(node1, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
sys1.name == Symbol(string(PlantModules.structmod(node1)) * string(PlantModules.id(node1)))

## get_MTK_system_dicts
MTK_system_dicts = PlantModules.get_MTK_system_dicts(graphs, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
length(MTK_system_dicts) == 2
length(MTK_system_dicts[1]) == 4
length(MTK_system_dicts[2]) == 1

## get_nb_nodes
node = node4
graphnr = 1
nb_nodes, nb_node_graphnrs = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
issetequal(nb_nodes, [node2, graphs[2]])

## get_connection_info
node = node4
graphnr = 1
PlantModules.get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, connecting_modules, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
#! test

## getnodeu0s
node = node4
structmodule = :Forest
func_module = lotka_volterra
nodeu0s = PlantModules.getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)
issetequal(nodeu0s, [:P => 0, :N => 5])


########################

sys = PlantModules.generate_system(model_defaults, module_defaults, module_coupling, struct_connections, func_connections,
	default_params = default_params, default_u0s = default_u0s
)

sys_simpl = structural_simplify(sys);

prob = ODEProblem(sys_simpl, [], (0, 48))
sol = solve(prob)

plot(sol)

PlantModules.plotnode(sol, PlantModules.nodes(graphs[2])[1])
PlantModules.plotnode(sol, PlantModules.nodes(graphs[1])[1], func_varname = :N)
PlantModules.plotgraph(sol, graphs[1])
PlantModules.plotgraph(sol, graphs[1], struct_module = :Forest)
PlantModules.plotgraph(sol, graphs[1], func_varname = :ΣF_P)
PlantModules.plotgraph(sol, graphs[1], struct_module = :Grassland, func_varname = :P)