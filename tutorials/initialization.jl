using Plots
using Pkg; Pkg.activate("./tutorials")
using ModelingToolkit, OrdinaryDiffEq

import ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel Foo begin
    @variables begin
        L(t)[1:3] = [1.0, 2.0, 3.0]
        V(t)
        P(t) = 1.0
    end

    @equations begin
        V ~ L[1]*L[2]*L[3]
        D(V) ~ P
        [D(L[i]) ~ P + D(P) for i in 1:3]...
    end
end

@mtkcompile sys = Foo()
prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob)
plot(sol)

using ForwardDiff

newprob = remake(prob, p = [sys.P => ForwardDiff.Dual(1.0, 1.0)])
newsol = solve(newprob)
plot(newsol)