# Structure #

Base.@kwdef mutable struct Forest <: PlantGraphs.Node
    N::Int = 20
    P::Int = 10
    δ::Float64 = 1.0
end

Base.@kwdef mutable struct Grassland <: PlantGraphs.Node
    N::Int = 100
    P::Int = 5
end

Base.@kwdef mutable struct Cave <: PlantGraphs.Node end

graph = Grassland() + Grassland(N = 50, P = 3) + (Forest(δ = 1.8), Forest(N = 5, P = 0, δ = 2.3))
graph2 = Cave() + Forest() + Grassland()
graphs = [graph, graph2]
intergraph_connections = [[1, 2] => (:Forest, :Cave)]
plantstructure = PlantStructure(graphs, intergraph_connections)

# Function #

## New functional modules

@independent_variables t
d = Differential(t);

function lotka_volterra(; name, α, β, δ, γ, N, P)
    @parameters α = α β = β δ = δ γ = γ
    @variables N(t) = N P(t) = P ΣF_N(t) ΣF_P(t)
    eqs = [
        d(N) ~ α*N - β*N*P + ΣF_N,
        d(P) ~ δ*N*P - γ*P + ΣF_P
    ]
    return System(eqs, t; name)
end

function fountain_of_rabbits(; name, η, N, P)
    @variables N(t) = N P(t) = P ΣF_N(t) ΣF_P(t)
    eqs = [
        d(N) ~ η + ΣF_N,
        d(P) ~ ΣF_P
    ]
    return System(eqs, t; name)
end

function wandering_animals(; name, κ)
    @parameters κ = κ
    @variables F_N(t) N_1(t) N_2(t) F_P(t) P_1(t) P_2(t)

    eqs = [
        F_N ~ κ * (N_2 - N_1),
        F_P ~ κ * (P_2 - P_1)
    ]

    wandering_eqs(node_MTK, nb_node_MTK, connection_MTK) = [
        connection_MTK.N_1 ~ node_MTK.N,
        connection_MTK.P_1 ~ node_MTK.P,
        connection_MTK.N_2 ~ nb_node_MTK.N,
        connection_MTK.P_2 ~ nb_node_MTK.P,
    ]

    return System(eqs, t; name), wandering_eqs
end

## Coupling

module_coupling = Dict(
    :Grassland => [lotka_volterra],
    :Forest => [lotka_volterra],
    :Cave => [fountain_of_rabbits]
)

connecting_modules = Dict(
    (:Grassland, :Forest) => wandering_animals,
    (:Forest, :Cave) => wandering_animals,
    (:Grassland, :Grassland) => wandering_animals
)

connecting_eqs(node_MTK, connection_MTKs) = [
    node_MTK.ΣF_N ~ sum([connection_MTK.F_N for connection_MTK in connection_MTKs]),
    node_MTK.ΣF_P ~ sum([connection_MTK.F_P for connection_MTK in connection_MTKs])
]

plantcoupling = PlantCoupling(; module_coupling, connecting_modules, connecting_eqs)

## Parameters

default_values = Dict(:α => 1.5, :β => 1.9, :δ => 1.5, :γ => 0.8, :N => 30, :P => 10, :κ => 0.01, :η => 10)

module_defaults = Dict(
    :Grassland => Dict(:β => 0.1),
    :Forest => Dict(:δ => 2.3),
    :Cave => Dict(:β => 0, :N => 1, :P => 0)
)

connection_values = Dict(
    (:Grassland, :Forest) => Dict(:κ => 0.03),
    (:Forest, :Cave) => Dict(:κ => 0.05),
    (:Grassland, :Grassland) => Dict(:κ => 0.2),
)

plantparams = PlantParameters(; default_values, module_defaults, connection_values)

# graphfuncs #

node1, node2, node3, node4 = collect(values(graph.nodes))

## getnodes
@test issetequal(PlantModules.getnodes(graph), [node1, node2, node3, node4])

## getneighbours
@test issetequal(PlantModules.getneighbours(node2, graph), [node1, node3, node4])

## getattributes
@test issetequal(PlantModules.getattributes(node4), [:P => 0, :N => 5, :δ => 2.3])

## getstructmod
@test PlantModules.getstructmod(node1) == :Grassland

## getid
@test allunique(PlantModules.getid.([node1, node2, node3, node4]))


# generate_system #

## getMTKsystem
node1, node2, node3, node4 = collect(values(graph.nodes))

checkunits = false
sys1 = PlantModules.getMTKsystem(node1, plantparams, plantcoupling, checkunits)
@test get_name(sys1) == Symbol(string(PlantModules.getstructmod(node1)) * string(PlantModules.getid(node1)))

## get_MTK_system_dicts
MTK_system_dicts = PlantModules.get_MTK_system_dicts(plantstructure, plantparams, plantcoupling, checkunits)
@test length(MTK_system_dicts) == 2
@test length(MTK_system_dicts[1]) == 4
@test length(MTK_system_dicts[2]) == 3

## get_nb_nodes
node = node4
graphnr = 1
nb_nodes, nb_node_graphnrs = PlantModules.get_nb_nodes(node, graphnr, plantstructure)
@test issetequal(nb_nodes, [node2, PlantModules.getnodes(graphs[2])[1]])

## get_connecting_module
node = node4
nb_node = nb_nodes[1]

connecting_module, reverse_order = PlantModules.get_connecting_module(node, nb_node, plantcoupling)

## get_connection_info
graphnr = 1
nb_node_graphnr = nb_node_graphnrs[1]

connection_MTK, connection_equations = PlantModules.get_connection_info(node, graphnr, nb_node, nb_node_graphnr, connecting_module,
	reverse_order, plantparams, MTK_system_dicts)

@test connection_MTK isa System
@test only(values(get_defaults(connection_MTK))) == 0.03
@test connection_equations isa Vector{Equation}

## getnodevalues
node = node4
structmodule = :Forest
func_module = lotka_volterra
nodevalues = PlantModules.getnodevalues(node, structmodule, func_module, plantparams)
@test issetequal(nodevalues, [:α => 1.5, :β => 1.9, :γ => 0.8,  :δ => 2.3, :N => 5, :P => 0])

node = PlantModules.getnodes(graph2)[1]
structmodule = :Cave
func_module = fountain_of_rabbits
nodeu0s = PlantModules.getnodevalues(node, structmodule, func_module, plantparams)
@test issetequal(nodeu0s, [:η => 10, :P => 0, :N => 1])

## generate_system
sys = PlantModules.generate_system(plantstructure, plantcoupling, plantparams)

prob = ODEProblem(sys, [], (0, 48))
sol = solve(prob)

# plot functions #
@test PlantModules.plotnode(sol, PlantModules.getnodes(graphs[2])[1]) isa Vector{AbstractPlot}
@test PlantModules.plotnode(sol, PlantModules.getnodes(graphs[1])[1], varname = :N) isa AbstractPlot

@test PlantModules.plotgraph(sol, graphs[1]) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, graphs[1], structmod = :Forest) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, graphs[1], varname = :ΣF_P) isa AbstractPlot
@test PlantModules.plotgraph(sol, graphs[1], structmod = :Grassland, varname = :P) isa AbstractPlot

@test PlantModules.plotgraph(sol, graphs) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, graphs, structmod = :Forest) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, graphs, varname = :ΣF_P) isa AbstractPlot
@test PlantModules.plotgraph(sol, graphs, structmod = :Grassland, varname = :P) isa AbstractPlot