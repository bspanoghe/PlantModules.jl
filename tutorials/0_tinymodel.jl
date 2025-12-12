using Revise, Infiltrator
using Plots
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using ModelingToolkit, OrdinaryDiffEq, PlantGraphs, Unitful
using Measurements

# testing

import PlantModules: t, d

function _hydraulic_module(; name, ϕ, E, Γ, V, Ψ)
    ρ_w = 1.0 # g / cm^3
    P = Ψ

    @constants (
        P_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / cm^3"],
    )
    @parameters (
        ϕ[1:2] = ϕ, [description = "Volumetric extensibility", unit = u"MPa^-1 * hr^-1"],
        E = E, [description = "Volumetric elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        W(t) = V * ρ_w, [description = "Water content", unit = u"g"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
    )

    eqs = [
        Ψ ~ P,
        W ~ V * ρ_w,
        
        d(P) ~ E * ((ΣF / ρ_w)/V - ϕ[1]*P_unit*logsumexp((P - Γ)/P_unit, α = 100)),
            # dV / V = ϕ * P + dP / E
        d(V) ~ ΣF / ρ_w,
    ]
    return System(eqs, t; name)
end

function _constant_hydraulic_connection(; name, K)
    @parameters (
        K = K, [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
    )

    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1)
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = [ 
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
    ]
    return System(eqs, t; name), get_connection_eqset
end


# structure

struct Compartment <: Node
    Ψ::Float64
    V::Float64
end

graph = Compartment(0.0, 10.0) + Compartment(-1.0, 3.0)
plantstructure = PlantStructure(graph)

# coupling
module_coupling = Dict(:Compartment => [_hydraulic_module])
connecting_modules = Dict((:Compartment, :Compartment) => _constant_hydraulic_connection)
plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# parameters
default_changes = Dict(:shape => Sphere(), :K => 1.0, :E => 1.0, :ϕ => [1.0, 0.0], :V => 0.0)
plantparams = PlantParameters(; default_changes)

# run it
tspan = (0.0, 10.0)
system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false);
prob = ODEProblem(system, [], tspan);
@time sol = solve(prob, saveat = 0.1);

plotgraph(sol, plantstructure, varname = :W)
plotgraph(sol, plantstructure, varname = :Ψ)

function re_solve(K)
    ps = get_subsystem_variables(system, plantstructure, :K, (:Compartment, :Compartment))
    newprob = remake(prob, p = Pair.(ps, [K]))

    return solve(newprob, saveat = 0.1)[system.Compartment1.V] .|> only
end

plot(re_solve(1.0))
plot(re_solve(1.1))


using ForwardDiff

K_grad = ForwardDiff.derivative(re_solve, 1.0)
plot(range(tspan..., step = 0.1), K_grad)