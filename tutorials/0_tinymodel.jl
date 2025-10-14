using Plots
using Revise, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using OrdinaryDiffEq, PlantGraphs

# for testing

using ModelingToolkit, Unitful
import ModelingToolkit: equations, parameters, observed
import PlantModules: t, d

function _hydraulic_module(; name, shape, ϕ_D, ϵ_D, Γ, T, D, Ψ, M, h)
    D, ϕ_D, ϵ_D = [PlantModules.correctdimensionality(shape, var) for var in [D, ϕ_D, ϵ_D]] 
        # turns scalar values into vectors of correct length

    num_D = getdimensionality(shape)
    
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
        ϕ_D[1:num_D] = ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        ϵ_D[1:num_D] = ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
        Pₕ = Pₕ, [description = "Gravitational water potential", unit = u"MPa"],
        g = g, [description = "Gravitational acceleration", unit = u"hN / g"] # (from N / kg) Pa = N/m^2 => MPa = hN/cm^2
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation ()
        W(t) = volume(shape, D) / ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        # V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"], 
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P + Π + Pₕ, # Water potential consists of a solute- and a pressure component
        Π ~ -R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        # V ~ W / ρ_w, # Volume is directly related to water content  
        # V ~ volume(shape, D), # Volume is also directly related to compartment dimensions
        W ~ ρ_w * volume(shape, D),
        [ΔP ~ ϵ_D[i] * (ΔD[i]/D[i] - ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100)) for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure
            # from ΔD[i]/D[i] ~ ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100) + ΔP/ϵ_D[i]

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return System(eqs, t; name)
end


# structure

struct Compartment <: Node
    Ψ::Float64
    D::Vector{Float64}
end

graph = Compartment(0.0, [10.0]) + Compartment(-1.0, [3.0])
plantstructure = PlantStructure(graph)

# coupling
module_coupling = Dict(:Compartment => [_hydraulic_module, constant_carbon_module, K_module])
connecting_modules = Dict((:Compartment, :Compartment) => hydraulic_connection)
plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# parameters
default_changes = Dict(:shape => Sphere(), :K => 1.0, :ϵ_D => 1.0)
plantparams = PlantParameters(; default_changes)

# run it

begin
    system = generate_system(plantstructure, plantcoupling, plantparams);
    prob = ODEProblem(system, [], (0.0, 10.0));
end;

isys = generate_initializesystem(system, checks = false)
isys_c = mtkcompile(isys, fully_determined = false)

equations(isys_c)
equations(isys)


@time sol = solve(prob);

plotgraph(sol, plantstructure, varname = :W)
plotgraph(sol, plantstructure, varname = :Ψ)