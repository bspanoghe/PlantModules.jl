# Function definition

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

# Coupling to structure

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

# Parameter definition

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

# generate_system
forest_node_1 = [node for node in getnodes(plantstructure) if getstructmod(node) == :Forest][1]
cave_node = [node for node in getnodes(plantstructure) if getstructmod(node) == :Cave] |> only

node = forest_node_1
nb_node = cave_node

## getMTKsystem

checkunits = false
sys1 = PlantModules.getMTKsystem(forest_node_1, plantparams, plantcoupling, checkunits)
@test get_name(sys1) == Symbol(string(PlantModules.getstructmod(forest_node_1)) * string(PlantModules.getid(forest_node_1)))

## get_MTK_system_dict
MTK_system_dict = PlantModules.get_MTK_system_dict(plantstructure, plantparams, plantcoupling, checkunits)
@test length(MTK_system_dict) == 7
@test MTK_system_dict[1] isa ModelingToolkit.System

## get_connecting_module
connecting_module, original_order = PlantModules.get_connecting_module(node, nb_node, plantcoupling)
@test connecting_module isa Function
@test original_order == true

## get_connection_info
connection_MTK, connection_equations = PlantModules.get_connection_info(node, nb_node, connecting_module,
	original_order, plantparams, MTK_system_dict)

@test connection_MTK isa ModelingToolkit.System
@test only(values(get_defaults(connection_MTK))) == 0.05
@test connection_equations isa Vector{Equation}

## getnodevalues
structmodule = :Forest
func_module = lotka_volterra
nodevalues = PlantModules.getnodevalues(node, structmodule, func_module, plantparams)
@test issetequal(nodevalues, [:α => 1.5, :β => 1.9, :γ => 0.8,  :δ => 1.8, :N => 20, :P => 10])

structmodule = :Cave
func_module = fountain_of_rabbits
nodeu0s = PlantModules.getnodevalues(nb_node, structmodule, func_module, plantparams)
@test issetequal(nodeu0s, [:η => 10, :P => 0, :N => 1])

## generate_system
sys = PlantModules.generate_system(plantstructure, plantcoupling, plantparams)

prob = ODEProblem(sys, [], (0, 48))
sol = solve(prob)