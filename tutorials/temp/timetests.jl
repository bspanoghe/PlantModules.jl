using Pkg; Pkg.activate("./tutorials")
using ModelingToolkit, DifferentialEquations, Plots

@variables t
d = Differential(t);

function hydraulic_compartment(; name, W)
    @variables W(t) = W ΔW(t) F_in(t) F_out(t)

    eqs = [
        ΔW ~ F_in - F_out,
        d(W) ~ ΔW,
    ]
    return ODESystem(eqs, t; name)
end

function hydraulic_connection(; name, K)
    @variables F(t) W_1(t) W_2(t)
    @parameters K = K

    eqs = [F ~ K*(W_2 - W_1)] # F is flux from compartment 2 to compartment 1

    return ODESystem(eqs, t; name)
end

function connect_compartments(compartment, nb, K)
    connection = hydraulic_connection(name = Symbol(string(compartment.name) * "nb"), K = K)
    
    connection_eqs = [
        connection.W_1 ~ compartment.W,
        connection.W_2 ~ nb.W,
        compartment.F_in ~ connection.F,
        compartment.F_out ~ 0
    ]

    return [connection], connection_eqs
end

function connect_compartments(compartment, nb_left, nb_right, K)
    connection_l = hydraulic_connection(name = Symbol(string(compartment.name) * "nb_left"), K = K)
    connection_r = hydraulic_connection(name = Symbol(string(compartment.name) * "nb_right"), K = K)

    connection_eqs = [
        connection_l.W_1 ~ compartment.W,
        connection_l.W_2 ~ nb_left.W,
        connection_r.W_1 ~ compartment.W,
        connection_r.W_2 ~ nb_right.W,
        compartment.F_in ~ connection_l.F,
        compartment.F_out ~ -connection_r.F
    ]

    return [connection_l, connection_r], connection_eqs
end

n_compartments = 300
K = 500.0
compartments = [hydraulic_compartment(name = Symbol("compartment$(i)"), W = i) for i in 1:n_compartments]
connections = ODESystem[]
connection_eqs = Equation[]

con, eqs = connect_compartments(compartments[1], compartments[2], K)
append!(connections, con)
append!(connection_eqs, eqs)

for comp_idx in 2:(n_compartments-1)
    con, eqs = connect_compartments(compartments[comp_idx], compartments[comp_idx-1], compartments[comp_idx+1], K)
    append!(connections, con)
    append!(connection_eqs, eqs)
end

con, eqs = connect_compartments(compartments[n_compartments], compartments[n_compartments - 1], K)
append!(connections, con)
append!(connection_eqs, eqs)

@named sys = ODESystem(connection_eqs, t, systems = vcat(compartments, connections))
sys_simpl = structural_simplify(sys);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 40))
@time sol = solve(prob);

# plot(sol)