using BenchmarkTools, Infiltrator, Revise
using Pkg; Pkg.activate("./tutorials")
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots

using ModelingToolkit: t_nounits as t, D_nounits as D

@connector function WaterPort(; W, name)
    pars = @parameters begin
        T = 293.15
        ρ = 1.0
    end
    vars = @variables begin
        W(t) = W
        F(t), [connect = Flow]
    end

    ODESystem(Equation[], t, vars, pars; name)
end

@component function System(; name)
    systems = @named begin
        w1 = WaterPort(; W = 100)
        w2 = WaterPort(; W = 0)
    end
    eqs = [
        connect(w1, w2)
        w1.F ~ w2.W - w1.W
        D(w1.W) ~ w1.F
        D(w2.W) ~ w2.F
    ]
    ODESystem(eqs, t, [], []; name, systems)
end

@named sys = System()
sys_simpl = structural_simplify(sys)
prob = ODEProblem(sys_simpl, [], (0, 10))
sol = solve(prob)
plot(sol)

#=
@mtkmodel WaterFall begin
    @components begin
        f = Water()
    end
    @equations begin
        f.W ~ 0
    end
end

@mtkmodel WaterPort begin
    @components begin
        p = Water()
        n = Water()
    end
    @variables begin
        ΔW(t)
        F(t)
    end
    @equations begin
        ΔW ~ p.W - n.W
        0 ~ p.F + n.F
        F ~ p.F
    end
end

@mtkmodel FlowingRes begin
    @extend WaterPort()
    @equations begin
        ΔW ~ F
    end
end

@mtkmodel FlowingCap begin
    @extend WaterPort()
    @equations begin
        D(ΔW) ~ F
    end
end

@mtkmodel Reservoir begin
    @extend WaterPort()
    @parameters begin
       Δw = 1.0
    end
    @equations begin
        Δw ~ ΔW
    end
end

@mtkmodel WaterCentral begin
    @components begin
        res = FlowingCap()
        cap = FlowingCap()
        source = Reservoir(Δw = 1.0)
        ground = WaterFall()
    end
    @equations begin
        connect(source.p, res.p)
        connect(res.n, cap.p)
        connect(cap.n, source.n)
        connect(cap.n, ground.f)
    end
end

@mtkbuild sys = WaterCentral()
prob = ODEProblem(sys, [sys.cap.ΔW => 0], (0, 10))
sol = solve(prob)
plot(sol)
plot(sol, idxs = [sys.res.p.W])
=#