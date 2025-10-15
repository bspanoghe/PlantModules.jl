using Plots
using Revise, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using OrdinaryDiffEq, PlantGraphs

# for testing

using ModelingToolkit, Unitful
import ModelingToolkit: equations, parameters, observed
import PlantModules: t, d

function _hydraulic_module(; name, D, Ψ, ϵ_D) 
    @parameters (
        ρ_w = 1.0, [description = "Density of water", unit = u"g / cm^3"],
        ϵ_D = ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
    )

    @variables (
        Ψ(t) = Ψ, [description = "Total water potential", unit = u"MPa"],
        W(t), [description = "Water content", unit = u"g"],
        D(t) = D, [description = "Radius", unit = u"cm"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔΨ(t), [description = "Change in water potential", unit = u"MPa / hr"],
        ΔD(t), [description = "Change in dimensions of compartment", unit = u"cm / hr", guess = 0.0],
    )

    eqs = [
        W ~ ρ_w * 4/3 * pi * D^3,
        ΔΨ ~ ϵ_D * ΔD/D,
        
        d(W) ~ ΣF,
        d(Ψ) ~ ΔΨ,
        d(D) ~ ΔD,
    ]
    return System(eqs, t; name, checks = false)
end


# structure

struct Compartment <: Node
    Ψ::Float64
    D::Float64
end

graph = Compartment(0.0, 10.0) + Compartment(-1.0, 3.0)
plantstructure = PlantStructure(graph)

# coupling
module_coupling = Dict(:Compartment => [_hydraulic_module])
connecting_modules = Dict((:Compartment, :Compartment) => constant_hydraulic_connection)
plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# parameters
default_changes = Dict(:shape => Sphere(), :K => 1000.0, :ϵ_D => 1.0)
plantparams = PlantParameters(; default_changes)

# run it

begin
    system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false);
    prob = ODEProblem(system, [], (0.0, 10.0), check_initialization_units = false);
end;

isys = generate_initializesystem(system, checks = false)
isys_c = mtkcompile(isys, fully_determined = false)

equations(isys_c)
equations(isys)


@time sol = solve(prob);

plotgraph(sol, plantstructure, varname = :W)
plotgraph(sol, plantstructure, varname = :Ψ)