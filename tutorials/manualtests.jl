using Revise #!
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Plots

function line(; name)
    @independent_variables t
    d = Differential(t)

    @variables (
        x(t) = 0.0,
        y(t)
    )
    eqs = [
        d(x) ~ y,
        y ~ cos(x)
    ]
    return ODESystem(eqs, t; name)
end

myline = line(name = :myline) |> structural_simplify
prob = ODEProblem(myline, [], (0.0, 10.0))
solve(prob) |> plot