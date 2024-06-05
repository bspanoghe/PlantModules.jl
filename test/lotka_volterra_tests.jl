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

default_params = (
    lotka_volterra = (α = 1.5, β = 1.9, δ = 1.5, γ = 0.8),
    fountain_of_rabbits = (η = 10,),
    wandering_animals = (κ = 0.01,)
)
default_u0s = (
    lotka_volterra = (N = 30, P = 10),
    fountain_of_rabbits = (N = 1, P = 0),
    wandering_animals = ()
)
module_defaults = (
    Grassland = (β = 0.1,),
    Forest = (δ = 2.3,),
    Cave = (β = 0,)
)
module_coupling = [
    lotka_volterra => [:Grassland, :Forest],
    fountain_of_rabbits => [:Cave]
]
graph = Grassland() + Grassland(N = 50, P = 3) + (Forest(δ = 1.8), Forest(N = 5, P = 0, δ = 2.3))
graph2 = Cave() + Forest() + Grassland()
graphs = [graph, graph2]
struct_connections = [graphs, [[1, 2] => (:Forest, :Cave)]] #! make vector of graphs as input unnecessary when theres only 1 graph
func_connections = [[() => wandering_animals], wandering_eqs]

graphs, intergraph_connections = struct_connections
connecting_modules, get_connection_eqs = func_connections


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
sys1 = PlantModules.getMTKsystem(node1, module_coupling, module_defaults, default_params, default_u0s, checkunits)
@test get_name(sys1) == Symbol(string(PlantModules.structmod(node1)) * string(PlantModules.id(node1)))

## get_MTK_system_dicts
MTK_system_dicts = PlantModules.get_MTK_system_dicts(graphs, module_coupling, module_defaults, default_params, default_u0s, checkunits)
@test length(MTK_system_dicts) == 2
@test length(MTK_system_dicts[1]) == 4
@test length(MTK_system_dicts[2]) == 3

## get_nb_nodes
node = node4
graphnr = 1
nb_nodes, nb_node_graphnrs = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
@test issetequal(nb_nodes, [node2, PlantModules.nodes(graphs[2])[1]])

## get_connection_info
node = node4
graphnr = 1
connection_MTKs, connection_equations = PlantModules.get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, connecting_modules, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
@test connection_MTKs isa Vector{ODESystem}
@test connection_equations isa Vector{Equation}

## getnodeparamu0s
node = node4
structmodule = :Forest
func_module = lotka_volterra
nodeu0s = PlantModules.getnodeparamu0s(node, structmodule, func_module, module_defaults, default_u0s)
@test issetequal(nodeu0s, [:P => 0, :N => 5])

node = PlantModules.nodes(graph2)[1]
structmodule = :Cave
func_module = fountain_of_rabbits
nodeu0s = PlantModules.getnodeparamu0s(node, structmodule, func_module, module_defaults, default_u0s)
@test issetequal(nodeu0s, [:P => 0, :N => 1])

## generate_system
sys = PlantModules.generate_system(default_params, default_u0s, module_defaults, module_coupling, struct_connections, func_connections)

sys_simpl = structural_simplify(sys);

prob = ODEProblem(sys_simpl, [], (0, 48))
sol = solve(prob)

# plot functions #
@test PlantModules.plotnode(sol, PlantModules.nodes(graphs[2])[1]) isa Vector
@test PlantModules.plotnode(sol, PlantModules.nodes(graphs[1])[1], func_varname = :N) isa Plots.Plot{T} where {T}

@test PlantModules.plotgraph(sol, graphs[1]) isa Vector
@test PlantModules.plotgraph(sol, graphs[1], struct_module = :Forest) isa Vector
@test PlantModules.plotgraph(sol, graphs[1], func_varname = :ΣF_P) isa Plots.Plot{T} where {T}
@test PlantModules.plotgraph(sol, graphs[1], struct_module = :Grassland, func_varname = :P) isa Plots.Plot{T} where {T}

@test PlantModules.plotgraph(sol, graphs) isa Vector
@test PlantModules.plotgraph(sol, graphs, struct_module = :Forest) isa Vector
@test PlantModules.plotgraph(sol, graphs, func_varname = :ΣF_P) isa Plots.Plot{T} where {T}
@test PlantModules.plotgraph(sol, graphs, struct_module = :Grassland, func_varname = :P) isa Plots.Plot{T} where {T}