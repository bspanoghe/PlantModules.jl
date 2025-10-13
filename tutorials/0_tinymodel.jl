using Plots
using Revise, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using OrdinaryDiffEq, PlantGraphs

# structure

struct Compartment <: Node
    Ψ::Float64
    D::Vector{Float64}
end

graph = Compartment(0.0, [10.0]) + Compartment(-1.0, [3.0])
plantstructure = PlantStructure(graph)

# coupling
module_coupling = Dict(:Compartment => [hydraulic_module, constant_carbon_module, K_module])
connecting_modules = Dict((:Compartment, :Compartment) => hydraulic_connection)
plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# parameters
default_changes = Dict(:shape => Sphere(), :K => 1.0, :ϵ_D => 1.0)
plantparams = PlantParameters(; default_changes)

# run it
system = generate_system(plantstructure, plantcoupling, plantparams)
prob = ODEProblem(system, [], (0.0, 10.0))
@time sol = solve(prob);

plotgraph(sol, plantstructure, varname = :W)
plotgraph(sol, plantstructure, varname = :Ψ)