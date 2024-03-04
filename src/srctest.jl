using Pkg; Pkg.activate(".")
include("PlantModules.jl")
using PlantGraphs, ModelingToolkit, Plots, DifferentialEquations, Unitful #! MTK imports etc. should not be necessary when package is done
import GLMakie.draw
using BenchmarkTools

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

function Ψ_soil_module(; name)
	@variables t Ψ(t) W_r(t)
	eqs = [Ψ ~ -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19)] #! check for realistic value of W_r
		# An empirical relationship between the soil water potential and relative water content

	return ODESystem(eqs; name)
end

function Ψ_air_module(; name) #! replace by constant Ψ?
	@variables t Ψ(t) W_r(t)
	@parameters T
	@constants R = 8.314e-6 #= MJ/K/mol =# V_w = 18e-6 #= m^3/mol =#

	eqs = [Ψ ~ R * T / V_w * log(W_r)]
		#! What's this equation called again? (Ask Jeroen)

	return ODESystem(eqs; name)
end

PlantModules.default_params
PlantModules.default_u0s

model_defaults = (Γ = 0.4, P = 0.2, T = 293.15) #! currently assumed variables with same name are the same for all functional modules

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
	Ψ_soil_module => [:Soil],
	Ψ_air_module => [:Air]
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
	:default => PlantModules.hydraulic_connection, #! make sure default works
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 8]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 4]), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-2]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-2]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

func_connections = [connecting_modules, get_connection_eqs]
#=
time_span = (0, 7*24.0) # We'll simulate our problem for a timespan of 7 days
prob = ODEProblem(plantsys, time_span)
sol = solve(prob)

plot(sol, struct_modules = [:Soil], func_vars = [:W]) #! imagine that this works
=#

include("PlantModules.jl")

system = PlantModules.get_system_definition(model_defaults, module_defaults, module_coupling,
	struct_connections, func_connections, checkunits = false)

	sys_simpl = structural_simplify(system)
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

mutable struct Foo <: Node
	bar::Int
end

mutable struct Oof <: Node
	rab::String
end

function qux(; name, bar, baz, quux, corge, xyzzy)
	@parameters bar = bar baz = baz quux = quux corge = corge
	@variables t xyzzy(t) = xyzzy
	eqs = [
		bar + baz ~ quux * xyzzy + corge
	]
	return ODESystem(eqs, t; name)
end

function snoo(; name, snee)
	@parameters snee = snee
	@variables t xyzzy1(t) xyzzy2(t) snaw(t)
	eqs = [
		snaw ~ xyzzy1 * xyzzy2 + snee
	]
	return ODESystem(eqs, t; name)
end

PlantModules.attributes(node::Foo) = Dict([:bar => node.bar])
PlantModules.nodetype(::Foo) = :Foo
PlantModules.nodetype(::Oof) = :Oof
PlantModules.id(::Foo) = 1

func_module = qux
structmodule = :Foo
func_connections = [:default => snoo]
module_coupling = [qux => [:Foo]]
default_params = (qux = (bar = 100, baz = 100, quux = 100, corge = 100), snoo = (snee = 100,))
default_u0s = (qux = (xyzzy = 100,),)
model_defaults = (baz = 42, quux = 42)
module_defaults = (Foo = (baz = 10, xyzzy = 10),)

## getnodeparams

node = Foo(1)

getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params) ==
	Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100])

## getnodeu0s

getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s) ==
	Dict([:xyzzy => 10])

## getMTKsystem

node = Foo(1)
testsystem = getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
testsystem.name == :(Foo1)
sort(Symbol.(keys(testsystem.defaults))) == sort(collect(keys(Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100, Symbol("xyzzy(t)") => 10]))))
sort(collect(values(testsystem.defaults))) == sort(collect(values(Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100, Symbol("xyzzy(t)") => 10]))))

## get_func_connection
node = Foo(1)
nb_node = Foo(2)
testsystem = get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
testsystem.name == :Foo1_Foo1
sort(Symbol.(keys(testsystem.defaults))) == sort(collect(keys(Dict([:snee => 100]))))
sort(collect(values(testsystem.defaults))) == sort(collect(values(Dict([:snee => 100]))))

## get_connection_info
graphnr = 1
graph = Foo(123) + (Foo(22), Foo(300))
node = PlantModules.nodes(graph)[1]
nb_nodes = PlantModules.nodes(graph)[2:3]
nb_node_graphnrs = [1, 1]
MTK_system_dicts = [Dict([PlantModules.nodes(graph)[nodenr].self_id =>
	qux(name = Symbol("Foo$nodenr"), bar = nodenr, baz = nodenr, quux = nodenr, corge = nodenr, xyzzy = nodenr)
	for nodenr in 1:3]
)]
get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
	[connection_MTK.xyzzy1 ~ node_MTK.xyzzy for connection_MTK in connection_MTKs]..., 
	[connection_MTK.xyzzy2 ~ nb_node_MTKs.xyzzy for (connection_MTK, nb_node_MTKs) in zip(connection_MTKs, nb_node_MTKs)]..., 
]

connection_MTKs, connection_equations = get_connection_info(node, graphnr, nb_nodes,
	nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts
)

connection_MTKs[1] == snoo(; name = Symbol("Foo" * string(PlantModules.id(node)) * "_" *
"Foo" * string(PlantModules.id(nb_nodes[1]))), snee = 100)
connection_MTKs[2] == snoo(; name = Symbol("Foo" * string(PlantModules.id(node)) * "_" *
"Foo" * string(PlantModules.id(nb_nodes[2]))), snee = 100)

connection_equations == get_connection_eqs(MTK_system_dicts[1][PlantModules.id(node)],
	[MTK_system_dicts[1][PlantModules.id(nb_nodes[1])], MTK_system_dicts[1][PlantModules.id(nb_nodes[2])]],
	connection_MTKs
)

## get_intergraph_neighbours
graphs = [Oof("1") + (Foo(2), Foo(3), Foo(4)), Oof("!")]
node_graphnr = 2
node = PlantModules.nodes(graphs[node_graphnr])[1]

### 1
intergraph_connections = [[1, 2] => (:Foo, :Oof)]
ig_neighbours, ig_neighbour_graphnrs = get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
ig_neighbour_graphnrs == [1, 1, 1]

### 2
intergraph_connections = [[2, 1] => (:Foo, :Oof)]
ig_neighbours, ig_neighbour_graphnrs = get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
isempty(ig_neighbours)
isempty(ig_neighbour_graphnrs)

### 3
intergraph_connections = [[2, 1] => (:Oof, :Foo)]
ig_neighbours, ig_neighbour_graphnrs = get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
ig_neighbour_graphnrs == [1, 1, 1]

## get_nb_nodes
graphs = [Oof("1") + (Foo(2), Foo(3), Foo(4)), Foo(90)]
graphnr = 1
node = PlantModules.nodes(graphs[graphnr])[1]
intergraph_connections = [[1, 2] => (:Oof, :Foo)]

nb_nodes = get_nb_nodes(node, graphs, graphnr, intergraph_connections)
nb_nodes[1] == vcat(PlantModules.nodes(graphs[1])[2:4], PlantModules.nodes(graphs[2])[1])
nb_nodes[2] == [1, 1, 1, 2]