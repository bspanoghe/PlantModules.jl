using Plots
using Revise, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using ModelingToolkit, OrdinaryDiffEq, PlantGraphs, Unitful
using Measurements

# testing

function _hydraulic_module(; name, ϕ, ϵ, Γ, T, V, Ψ, M, h)
    R = 8.314 # MPa * cm^3 / K / mol
    ρ_w = 1.0 # g / cm^3
    g = 9.8 * 1e-5 # hN / g
    Pₕ = ρ_w * g * h # MPa
    P = Ψ + R*T*M - Pₕ # MPa

    @constants (
        R = R, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
        P_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / cm^3"],
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϕ = ϕ, [description = "Volumetric extensibility", unit = u"MPa^-1 * hr^-1"],
        ϵ = ϵ, [description = "Volumetric elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
        Pₕ = Pₕ, [description = "Gravitational water potential", unit = u"MPa"],
        g = g, [description = "Gravitational acceleration", unit = u"hN / g"] # (from N / kg) Pa = N/m^2 => MPa = hN/cm^2
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation ()
        W(t) = ρ_w * V, [description = "Water content", unit = u"g"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
    )

    eqs = [
        Ψ ~ P + Π + Pₕ,
        Π ~ -R*T*M,
        W ~ ρ_w * V,
        
        d(P) ~ ϵ * (ΣF/V - ϕ*P_unit*logsumexp((P - Γ)/P_unit, α = 100)),
        d(V) ~ ΣF,
    ]
    return System(eqs, t; name)
end



# structure

struct Compartment <: Node
    Ψ::Float64
    V::Float64
end

graph = Compartment(0.0, 10.0) + Compartment(-1.0, 3.0)
plantstructure = PlantStructure(graph)

# coupling
module_coupling = Dict(:Compartment => [_hydraulic_module, constant_carbon_module])
connecting_modules = Dict((:Compartment, :Compartment) => constant_hydraulic_connection)
plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# parameters
default_changes = Dict(:shape => Sphere(), :K => 1000.0, :ϵ => 1.0, :ϕ => 1.0)
plantparams = PlantParameters(; default_changes)

# run it
system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false);
prob = ODEProblem(system, Float64[], (0.0, 10.0), check_initialization_units = false);
@time sol = solve(prob, saveat = 0.1);

plotgraph(sol, plantstructure, varname = :W)
plotgraph(sol, plantstructure, varname = :D)

plotgraph(sol, plantstructure, varname = :Ψ)

function resolve(K)
    ps = get_subsystem_variables(system, plantstructure, :K, (:Compartment, :Compartment))
    newprob = remake(prob, p = Pair.(ps, [K]))
    return solve(newprob, saveat = 0.1)[system.Compartment1.D] .|> only
end

plot(resolve(10.0))
plot(resolve(1.0))

using ForwardDiff

ForwardDiff.derivative(resolve, 1.0)