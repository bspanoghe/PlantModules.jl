using BenchmarkTools, Infiltrator, Revise
using Pkg; Pkg.activate("./tutorials")
using ModelingToolkit, OrdinaryDiffEq, Unitful
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

function WaterPort(; name, W)
    pars = @parameters begin
        T = 298.15
    end
    vars = @variables begin
        W(t) = W
        F(t)#, [connect = Flow]
    end
    eqs = [
        D(W) ~ F
    ]

    ODESystem(eqs, t, vars, pars; name)
end

function WaterLoo(; name)


end

function System(; name)
    systems = @named begin
        w1 = WaterPort(; W = 100)
        w2 = WaterPort(; W = 0)
        w3 = WaterPort(; W = 50)
    end
    eqs = [
        w1.F ~ w2.W - w1.W
        w2.F ~ -w1.F
        # D(w1.W) ~ w1.F
        # D(w2.W) ~ w2.F
    ]
    ODESystem(eqs, t, [], []; name, systems)
end

@named sys = System()
sys_simpl = structural_simplify(sys)
prob = ODEProblem(sys_simpl, [], (0, 10))
sol = solve(prob)
plot(sol)











@connector Water begin
    W(t) = W
    F(t), [connect = Flow]
end

@mtkmodel Connection begin
    @components begin
        p = Water()
        n = Water()
    end
    @variables begin
        F(t)
    end
    @equations begin
        F ~ p.W - n.W
        p.F ~ -n.F
        F ~ p.F
    end
end

@component function System(; name) begin
    @components begin
        w1_2 = Water(W = 100)
        w1_3 = Water(W = 100)
        w2 = Water(W = 0)
        w3 = Water(W = 50)
        con1 = Connection()
        con2 = Connection()
    end

    @equations begin
        connect(w1_2, w1_3)

        connect(w1_2, con1.p)
        connect(w2, con1.n)

        connect(w1_3, con2.p)
        connect(w3, con2.n)

        w1_2.W ~ 100
    end
end

@mtkbuild sys = System()
prob = ODEProblem(sys, [], (0, 10))
sol = solve(prob)
plot(sol)



function water_comp(; name, W)
    pars = @parameters (
        T = 293.15
    )
    vars = @variables (
        W(t) = W,
        ΣF(t)
    )
    eqs = [
        D(W) ~ ΣF
    ]

    ODESystem(eqs, t, vars, pars; name)
end

function connection(; name, K, systems)
    pars = @parameters (
        K = K
    )
    vars = @variables (
        F(t)
    )
    eqs = [
        D(W) ~ ΣF
    ]

    ODESystem(eqs, t, vars, pars; name)
end