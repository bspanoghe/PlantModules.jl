using BenchmarkTools, Infiltrator, Revise
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots

using ModelingToolkit: t_nounits as t, D_nounits as D

@connector Water begin
    F(t), [connect = Flow]
    W(t)
end

@mtkmodel WComp begin
    @extend Water()
    @variables begin
        Ψ(t)
        W(t)
    end
    @equations begin
        D(W) ~ F
    end
end

@named hoo = WComp(Ψ = 3)

@mtkmodel WaterPort begin
    @components begin
        c1 = WComp()
        c2 = WComp()
    end
    @equations begin
        0 ~ c1.F + c2.F
        c2.Ψ - c1.Ψ ~ c1.F
    end
end

@mtkmodel Growing begin
    @components begin
        w = Water()
    end
    @variables begin
        V(t)        
    end
    @equations begin
        
    end
end

@mtkmodel Env begin
    @components begin
        w = Water()
    end
    @equations begin
        w.Ψ ~ 0
    end
end

@named sys = System()
sys_simpl = structural_simplify(sys)
prob = ODEProblem(sys_simpl, [], (0, 10))
sol = solve(prob)

plot(sol)