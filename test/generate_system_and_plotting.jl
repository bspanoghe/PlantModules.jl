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
struct_connections = PlantStructure(graphs, intergraph_connections)

# Function #

## New functional modules

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

    wandering_eqs(node_MTK, nb_node_MTK, connection_MTK, reverse_order) = [
        connection_MTK.N_1 ~ node_MTK.N,
        connection_MTK.P_1 ~ node_MTK.P,
        connection_MTK.N_2 ~ nb_node_MTK.N,
        connection_MTK.P_2 ~ nb_node_MTK.P,
    ]

    return ODESystem(eqs, t; name), wandering_eqs
end

default_values = Dict(
    lotka_volterra => Dict(:α => 1.5, :β => 1.9, :δ => 1.5, :γ => 0.8, :N => 30, :P => 10),
    fountain_of_rabbits => Dict(:η => 10, :N => 1, :P => 0),
    wandering_animals => Dict(:κ => 0.01)
)

module_defaults = Dict(
    :Grassland => Dict(:β => 0.1),
    :Forest => Dict(:δ => 2.3),
    :Cave => Dict(:β => 0)
)

connecting_modules = [
    (:Grassland, :Forest) => (wandering_animals, Dict(:κ => 0.03)),
    (:Forest, :Cave) => (wandering_animals, Dict(:κ => 0.05)),
    (:Grassland, :Grassland) => (wandering_animals, Dict(:κ => 0.2)),
]

multi_connection_eqs(node_MTK, connection_MTKs) = [
    node_MTK.ΣF_N ~ sum([connection_MTK.F_N for connection_MTK in connection_MTKs]),
    node_MTK.ΣF_P ~ sum([connection_MTK.F_P for connection_MTK in connection_MTKs])
]

func_connections = PlantFunctionality(default_values = default_values, module_defaults = module_defaults,
    connecting_modules = connecting_modules, connecting_eqs = multi_connection_eqs)

module_coupling = Dict(
    :Grassland => [lotka_volterra],
    :Forest => [lotka_volterra],
    :Cave => [fountain_of_rabbits]
)

# graphfuncs #

node1, node2, node3, node4 = collect(values(graph.nodes))

## nodes
@test issetequal(PlantModules.nodes(graph), [node1, node2, node3, node4])

## neighbours
@test issetequal(PlantModules.neighbours(node2, graph), [node1, node3, node4])

## attributes
@test issetequal(PlantModules.attributes(node4), [:P => 0, :N => 5, :δ => 2.3])

## structmod
@test PlantModules.structmod(node1) == :Grassland

## id
@test allunique(PlantModules.id.([node1, node2, node3, node4]))


# generate_system #

## getMTKsystem
node1, node2, node3, node4 = collect(values(graph.nodes))

checkunits = false
sys1 = PlantModules.getMTKsystem(node1, default_values, module_defaults, module_coupling, checkunits)
@test get_name(sys1) == Symbol(string(PlantModules.structmod(node1)) * string(PlantModules.id(node1)))

## get_MTK_system_dicts
MTK_system_dicts = PlantModules.get_MTK_system_dicts(graphs, default_values, module_defaults, module_coupling, checkunits)
@test length(MTK_system_dicts) == 2
@test length(MTK_system_dicts[1]) == 4
@test length(MTK_system_dicts[2]) == 3

## get_nb_nodes
node = node4
graphnr = 1
nb_nodes, nb_node_graphnrs = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
@test issetequal(nb_nodes, [node2, PlantModules.nodes(graphs[2])[1]])

## get_connecting_module
node = node4
nb_node = nb_nodes[1]

connecting_module, reverse_order = PlantModules.get_connecting_module(node, nb_node, connecting_modules)

## get_connection_info
graphnr = 1
nb_node_graphnr = nb_node_graphnrs[1]

connection_MTK, connection_equations = PlantModules.get_connection_info(node, graphnr, nb_node, nb_node_graphnr,
 connecting_module, reverse_order, default_values, MTK_system_dicts
)

@test connection_MTK isa ODESystem
@test only(values(connection_MTK.defaults)) == 0.03
@test connection_equations isa Vector{Equation}

## getnodevalues
node = node4
structmodule = :Forest
func_module = lotka_volterra
nodevalues = PlantModules.getnodevalues(node, structmodule, func_module, module_defaults, default_values)
@test issetequal(nodevalues, [:α => 1.5, :β => 0.1, :γ => 0.8,  :δ => 2.3, :N => 5, :P => 0])

node = PlantModules.nodes(graph2)[1]
structmodule = :Cave
func_module = fountain_of_rabbits
nodeu0s = PlantModules.getnodevalues(node, structmodule, func_module, module_defaults, default_values)
@test issetequal(nodeu0s, [:η => 10, :P => 0, :N => 1])

## generate_system
sys = PlantModules.generate_system(struct_connections, func_connections, module_coupling)

sys_simpl = structural_simplify(sys);

prob = ODEProblem(sys_simpl, [], (0, 48))
sol = solve(prob)

# plot functions #
@test PlantModules.plotnode(sol, PlantModules.nodes(graphs[2])[1]) isa Vector
@test PlantModules.plotnode(sol, PlantModules.nodes(graphs[1])[1], varname = :N) isa Plots.Plot{T} where {T}

@test PlantModules.plotgraph(sol, graphs[1]) isa Vector
@test PlantModules.plotgraph(sol, graphs[1], structmod = :Forest) isa Vector
@test PlantModules.plotgraph(sol, graphs[1], varname = :ΣF_P) isa Plots.Plot{T} where {T}
@test PlantModules.plotgraph(sol, graphs[1], structmod = :Grassland, varname = :P) isa Plots.Plot{T} where {T}

@test PlantModules.plotgraph(sol, graphs) isa Vector
@test PlantModules.plotgraph(sol, graphs, structmod = :Forest) isa Vector
@test PlantModules.plotgraph(sol, graphs, varname = :ΣF_P) isa Plots.Plot{T} where {T}
@test PlantModules.plotgraph(sol, graphs, structmod = :Grassland, varname = :P) isa Plots.Plot{T} where {T}