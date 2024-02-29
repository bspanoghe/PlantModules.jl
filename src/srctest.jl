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

C_root = 800 # (mol/m^3) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 1200 # Idem here
C_leaf = 1800 # And here!
# loosely based on https://www.researchgate.net/figure/Soluble-sugar-content-A-starch-content-B-non-structural-carbon-NSC-content-C_fig1_362956616

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
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.7, 0.1, 0.05]), D = [0.3, 0.05, 0.03], M = C_root), #! changed C to M_const #! include check whether module_defaults variabkle names correspond with default_params / default_u0s ?
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.8, 0.03]), D = [0.4, 0.03],  M = C_stem),
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

plant_graph = Root() + Stem() + (Leaf([0.25]), Leaf([0.15]))

# draw(plant_graph, resolution = (500, 400))

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = [graphs, intergraph_connections]

func_connections = [
	:default => PlantModules.hydraulic_connection, #! make sure default works
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 50]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 800]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-3]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-2]) #! check value
]

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

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

###################### source code begin ######################

# Get MTK system corresponding with node
function getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
	structmodule = PlantModules.nodetype(node)
	func_modules = [coupling.first for coupling in module_coupling if structmodule in coupling.second]

	component_systems = Vector{ODESystem}(undef, length(func_modules))

	for (modulenum, func_module) in enumerate(func_modules)
		nodeparams = getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params)
		nodeu0s = getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)

		component_systems[modulenum] = func_module(; :name => :foo, # real name given later
			Pair.(keys(nodeparams), values(nodeparams))..., Pair.(keys(nodeu0s), values(nodeu0s))...)
	end

	MTKsystem = collapse(component_systems, name = Symbol(string(structmodule) * string(PlantModules.id(node))))

	return MTKsystem
end

# get correct parameter values for node between those defined in the model defaults, module defaults and node values
function getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params)
	if !haskey(default_params, Symbol(func_module))
		return Dict() # no parameters specified in default_params => functional module has no parameters
	end

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
	if !haskey(default_u0s, Symbol(func_module))
		return Dict() # no u0s specified in default_u0s => functional module has no u0s
	end
	
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

# collapse multiple ODESystems into one. like ModelingToolkit.compose, but keeps a single namespace
function collapse(systems::Vector{ODESystem}; name::Symbol)
    return ODESystem(vcat([system.eqs for system in systems]...), systems[1].iv, vcat([states(system) for system in systems]...),
        vcat([parameters(system) for system in systems]...), name = name)
end

# get MTK systems and the vector of equations defining connections between MTK systems of a node and its neighbour nodes
function get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
	node_MTK = MTK_system_dicts[graphnr][PlantModules.id(node)]

	nb_node_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
	connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes))

	for (nb_nr, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs))
		nb_node_MTKs[nb_nr] = MTK_system_dicts[nb_node_graphnr][PlantModules.id(nb_node)]
		connection_MTKs[nb_nr] = get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
	end

	connection_equations = get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs)
	return connection_MTKs, connection_equations
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
	
	node_intergraph_connections_plus_igidxs = [(intergraph_connection, findfirst(node_graphnr .== intergraph_connection[1]))
		for intergraph_connection in intergraph_connections 
		if node_graphnr in intergraph_connection[1]
	] # igidx is index of the node's graph in the intergraph connection

	for (node_intergraph_connection, igidx) in node_intergraph_connections_plus_igidxs
		if PlantModules.nodetype(node) == node_intergraph_connection[2][igidx]
			nb_igidx = 3-igidx # 2 if igidx == 1 and 1 if igidx == 2
			nb_graphnr = node_intergraph_connection[1][nb_igidx]
			nb_graph = graphs[nb_graphnr]
			connection_nb_nodes = [nb_node for nb_node in PlantModules.nodes(nb_graph) if PlantModules.nodetype(nb_node) == node_intergraph_connection[2][nb_igidx]]
			append!(nb_nodes, connection_nb_nodes)
			append!(nb_node_graphnrs, repeat([nb_graphnr], length(connection_nb_nodes)))
		end
	end

	return nb_nodes, nb_node_graphnrs
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

# get node idx to corresponding MTK system
get_MTK_system_dicts(graphs, module_coupling, module_defaults, model_defaults, default_params, default_u0s) = 
	[
		[PlantModules.id(node) => 
			getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
			for node in PlantModules.nodes(graph)
		] |> Dict for graph in graphs
]

###################### source code end ######################



###################### main function begin ######################

function get_system_definition(model_defaults, module_defaults,	module_coupling,
	struct_connections, func_connections; checkunits = true
	)

	default_params = PlantModules.default_params # is this needed here?
	default_u0s = PlantModules.default_u0s
	
	graphs = struct_connections[1]

	MTK_system_dicts = get_MTK_system_dicts(graphs, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
	MTK_systems = vcat(collect.(values.(MTK_system_dicts))...)

	connection_MTKs = ODESystem[]
	connection_equations = Equation[]

	for (graphnr, graph) in enumerate(graphs)
		for node in PlantModules.nodes(graph)
			nb_nodes, nb_node_graphnrs = get_nb_nodes(node, graphs, graphnr, intergraph_connections)

			if !isempty(nb_nodes) # no neighbours no connection info
				node_connection_MTK, node_connection_equations = get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
				append!(connection_MTKs, node_connection_MTK)
				append!(connection_equations, node_connection_equations)
			end
		end
	end

	system = ODESystem(connection_equations, name = :system, systems = vcat(MTK_systems, connection_MTKs), checks = checkunits)

	return system
end
###################### main function end ######################

system = get_system_definition(model_defaults, module_defaults,	module_coupling,
	struct_connections, func_connections, checkunits = false)
sys_simpl = structural_simplify(system)
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 48))
sol = solve(prob)

plot(sol, idxs = [Symbol("Soil1₊Ψ"), Symbol("Root1₊Ψ"), Symbol("Stem2₊Ψ"), Symbol("Leaf3₊Ψ"), Symbol("Air1₊Ψ")])
plot(sol, idxs = [Symbol("Root1₊Π"), Symbol("Stem2₊Π"), Symbol("Leaf3₊Π")])
plot(sol, idxs = [Symbol("Root1₊P"), Symbol("Stem2₊P"), Symbol("Leaf3₊P")])

plot(sol, idxs = [Symbol("Root1₊M_amount")])


plot(sol, idxs = [Symbol("Soil1₊W")])
plot(sol, idxs = [Symbol("Root1₊W")])
plot(sol, idxs = [Symbol("Stem2₊W")])
plot(sol, idxs = [Symbol("Leaf3₊W")])
plot(sol, idxs = [Symbol("Air1₊W")])


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