using Pkg; Pkg.activate(".")
include("PlantModules.jl")
using PlantGraphs, ModelingToolkit, GLMakie, DifferentialEquations #! MTK imports etc. should not be necessary when package is done
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


C_root = 0.5 # We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 1 # Idem here
C_leaf = 3 # And here!
Ψ_soil_func(W_r) = -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19) # An empirical relationship between the soil water potential and relative water content
Ψ_air_func(T, W_r) = R * T / V_w * log(W_r) #! What's this equation called again? (Ask Jeroen)

@constants R = 8.314 V_w = 18e-6

PlantModules.default_params
PlantModules.default_u0s

model_defaults = (Γ = 0.4, P = 0.2, M = 15.0, T = 293.15) #! currently assumed variables with same name are the same for all functional modules

module_defaults = (
	Root = (D = [0.3, 0.05, 0.03], shape = PlantModules.Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.7, 0.1, 0.05]), M = C_root), #! changed C to M #! include check whether module_defaults variabkle names correspond with default_params / default_u0s ?
	Stem = (D = [0.4, 0.03], shape = PlantModules.Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.8, 0.03]), M = C_stem),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 500.0, T = 288.15, Ψ = Ψ_soil_func),
	Air = (Ψ = Ψ_air_func,)
)

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air]
]

plant_graph = Root() + Stem() + (Leaf([0.25]), Leaf([0.15]))

# draw(plant_graph, resolution = (500, 400))

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]


intergraph_connections = [[1, 2] => (:Soil, :Root), [1, 3] => (:Air, :Leaf)]

graphs = [Stem() + Leaf([0.3])] #! remove
intergraph_connections = [] #! remove

struct_connections = [graphs, intergraph_connections]

func_connections = [
	:default => PlantModules.hydraulic_connection,
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 50]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 800]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-3]) 
]

plantsys = PlantModules.PlantSystem(
	model_defaults = model_defaults,
	module_defaults = module_defaults,
	module_coupling = module_coupling,
	struct_connections = struct_connections,
	func_connections = func_connections
)

#=
time_span = (0, 7*24.0) # We'll simulate our problem for a timespan of 7 days
prob = ODEProblem(plantsys, time_span)
sol = solve(prob)

plot(sol, struct_modules = [:Soil], func_vars = [:W]) #! imagine that this works
=#

#######

# Get MTK system corresponding with node
function getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
	structmodule = PlantModules.nodetype(node)
	func_module = [coupling.first for coupling in module_coupling if structmodule in coupling.second][1]

	nodeparams = getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params)
	nodeu0s = getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)

	MTKsystem = func_module(; :name => Symbol(string(structmodule) * string(PlantModules.id(node))),
		Pair.(keys(nodeparams), values(nodeparams))..., Pair.(keys(nodeu0s), values(nodeu0s))...)
	return MTKsystem
end

# get correct parameter values for node between those defined in the model defaults, module defaults and node values
function getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params)
	params = default_params[Symbol(func_module)] |> x -> Dict(Pair.(keys(x), values(x))) # PlantModules defaults
	paramnames = keys(params)

	node_attributes = PlantModules.attributes(node)
	nodemoduledefaults = module_defaults[structmodule]

	for paramname in paramnames
		if paramname in keys(node_attributes) # Change to node-specific value
			params[paramname] = node_attributes[paramname]
		elseif paramname in keys(nodemoduledefaults) # Change to module-wide defaults
			params[paramname] = nodemoduledefaults[paramname]
		elseif paramname in keys(model_defaults) # Change to model-wide defaults
			params[paramname] = model_defaults[paramname]
		end
	end

	return params
end

# get correct inital variable values for node between those defined in the model defaults, module defaults and node values
function getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)
	u0s = default_u0s[Symbol(func_module)] |> x -> Dict(Pair.(keys(x), values(x))) # PlantModules defaults
	u0names = keys(u0s)

	node_attributes = PlantModules.attributes(node)
	nodemoduledefaults = module_defaults[structmodule]

	for u0name in u0names
		if u0name in keys(node_attributes) # Change to node-specific value
			u0s[u0name] = node_attributes[u0name]
		elseif u0name in keys(nodemoduledefaults) # Change to module-wide defaults
			u0s[u0name] = nodemoduledefaults[u0name]
		elseif u0name in keys(model_defaults) # Change to model-wide defaults
			u0s[u0name] = model_defaults[u0name]
		end
	end

	return u0s
end

# get MTK systems and the vector of equations defining connections between MTK systems of a node and its neighbour nodes
function get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
	node_MTK = MTK_system_dicts[graphnr][PlantModules.id(node)]

	nb_node_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
	connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
	connection_keys = Vector{Tuple{Tuple, Tuple}}(undef, length(nb_nodes))

	for (nb_nr, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs))
		nb_node_MTKs[nb_nr] = MTK_system_dicts[nb_node_graphnr][PlantModules.id(nb_node)]
		connection_MTKs[nb_nr] = get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
		connection_keys[nb_nr] = ((graphnr, PlantModules.id(node)), (nb_node_graphnr, PlantModules.id(nb_node)))
	end

	connection_equations = get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs)
	connection_dict = Dict(Pair.(connection_keys, connection_MTKs))
	return connection_dict, connection_equations
end

# given two nodes' structural modules, get the MTK system of the functional connection between them
function get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
	structmodule = PlantModules.nodetype(node)
	nb_structmodule = PlantModules.nodetype(nb_node)

	func_connection_keys = first.(func_connections)
	connection_idx = findfirst([(structmodule, nb_structmodule)] .== func_connection_keys .|| [(nb_structmodule, structmodule)] .== func_connection_keys) # assuming symmetry in functional connections
	
	if isnothing(connection_idx)
		connector_func, connection_specific_values = func_connections[1][2], []
	else
		connector_func, connection_specific_values = func_connections[connection_idx].second
	end

	default_conn_info = merge(get(default_params, Symbol(connector_func), []), get(default_u0s, Symbol(connector_func), []))
	conn_info = merge(default_conn_info, connection_specific_values)

	func_connection = connector_func(;
		name = Symbol(string(structmodule) * string(PlantModules.id(node)) * "_" *
			string(nb_structmodule) * string(PlantModules.id(nb_node))),
		Pair.(keys(conn_info), values(conn_info))...
	)

	return func_connection
end

function get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
	nb_nodes = []
	nb_node_graphnrs = []
	
	node_intergraph_connections = [intergraph_connection for intergraph_connection in intergraph_connections if node_graphnr in intergraph_connection[1]]
	for node_intergraph_connection in node_intergraph_connections
		if PlantModules.nodetype(node) in node_intergraph_connection[2]
			nb_graphnr = node_intergraph_connection[1][node_intergraph_connection[1] .!= node_graphnr][1]
			nb_graph = graphs[nb_graphnr]
			connection_nb_nodes = [nb_node for nb_node in PlantModules.nodes(nb_graph) if PlantModules.nodetype(nb_node) in node_intergraph_connection[2]]
			append!(nb_nodes, connection_nb_nodes)
			append!(nb_node_graphnrs, repeat([nb_graphnr], length(connection_nb_nodes)))
		end
	end

	return nb_nodes, nb_node_graphnrs # also return graph nrs for looking up MTK things?
end

# get neighbouring nodes of a node both from the same graph and all connected graphs
function get_nb_nodes(node, graphs, graphnr, intergraph_connections)
	graph = graphs[graphnr]
	
	intra_nb_nodes = PlantModules.neighbours(node, graph)
	intra_nb_node_graphnrs = repeat([graphnr], length(intra_nb_nodes))

	inter_nb_nodes, inter_nb_node_graphnrs = get_intergraph_neighbours(node, graphnr, graphs, intergraph_connections)

	nb_nodes = vcat(intra_nb_nodes, inter_nb_nodes)
	nb_node_graphnrs = vcat(intra_nb_node_graphnrs, inter_nb_node_graphnrs)

	return nb_nodes, nb_node_graphnrs
end

function get_symmetry_info(connection_dict, get_symmetry_eqs)
	symmetry_eqs = Equation[]
	for connection_key in keys(connection_dict)
		rev_key = reverse(connection_key) # key of reverse connection (by construction of keys)
		append!(symmetry_eqs, get_symmetry_eqs(connection_dict[connection_key], connection_dict[rev_key]))
	end

	return symmetry_eqs
end

###################### source code end ######################

graphs = struct_connections[1]
graphnr = 1
graph = graphs[graphnr]
node = PlantModules.nodes(graph)[1]
default_params = PlantModules.default_params
default_u0s = PlantModules.default_u0s

get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
	[connection_MTK.Ψ_1 ~ node_MTK.Ψ for connection_MTK in connection_MTKs]...,
	[connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ for (connection_MTK, nb_node_MTK) in zip(connection_MTKs, nb_node_MTKs)]...,
	node_MTK.ΣF ~ sum([connection_MTK.F for connection_MTK in connection_MTKs])
] #! put this in func_modules and explain to user they gotta provide this stuff mr white yo

get_symmetry_eqs(connection_MTK, reverse_connection_MTK) = [
	connection_MTK.F ~ -reverse_connection_MTK.F
] #! put this in func_modules and explain to user they gotta provide this stuff mr white yo




###################### main function begin ######################
connection_dict = Dict{Tuple{Tuple, Tuple}, ODESystem}()
connection_equations = Equation[]

MTK_system_dicts = [[PlantModules.id(node) => getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
	for node in PlantModules.nodes(graph)] |> Dict for graph in graphs] # node idx to corresponding MTK system

MTK_systems = vcat(collect.(values.(MTK_system_dicts))...)

for (graphnr, graph) in enumerate(graphs)
	for node in PlantModules.nodes(graph)

		nb_nodes, nb_node_graphnrs = get_nb_nodes(node, graphs, graphnr, intergraph_connections)

		if !isempty(nb_nodes) # no neighbours no connection info
			node_connection_dict, node_connection_equations = get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
			merge!(connection_dict, node_connection_dict)
			append!(connection_equations, node_connection_equations)
		end

	end
end

append!(connection_equations, get_symmetry_info(connection_dict, get_symmetry_eqs))
connection_MTKs = collect(values(connection_dict))
###################### main function end ######################

system = ODESystem(connection_equations, name = :system, systems = vcat(MTK_systems, connection_MTKs)) |> structural_simplify
prob = ODEProblem(system, ModelingToolkit.missing_variable_defaults(system), (0.0, 10))
sol = solve(prob) # i dont think its valid

"""
ArgumentError: SymbolicUtils.BasicSymbolic{Real}[var"Root5₊D(t)[1]ˍt", var"Stem6₊D(t)[1]ˍt", var"Leaf7₊D(t)[1]ˍt",
var"Leaf8₊D(t)[1]ˍt", Soil1_Root5₊Ψ_2(t), Air1_Leaf8₊Ψ_2(t), Leaf7_Stem6₊Fˍt(t), Soil1_Root5₊Ψ_2ˍt(t), Air1_Leaf8₊Fˍt(t),
Air1_Leaf7₊Fˍt(t), Leaf8_Stem6₊Ψ_1ˍt(t), Leaf8_Stem6₊Fˍt(t)] are missing from the variable map.
"""

# tests #

mutable struct Foo <: Node
	bar::Int
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
graph = Foo(1) + (Foo(2), Foo(3))
node = PlantModules.nodes(graph)[1]
MTK_system_dict = Dict([PlantModules.nodes(graph)[nodenr].self_id =>
	qux(name = Symbol("Foo$nodenr"), bar = nodenr, baz = nodenr, quux = nodenr, corge = nodenr, xyzzy = nodenr)
	for nodenr in 1:3]
)
get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
	[connection_MTK.xyzzy1 ~ node_MTK.xyzzy for connection_MTK in connection_MTKs]..., 
	[connection_MTK.xyzzy2 ~ nb_node_MTKs.xyzzy for (connection_MTK, nb_node_MTKs) in zip(connection_MTKs, nb_node_MTKs)]..., 
]

connection_MTKs, connection_equations = get_connection_info(node, graph, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dict)