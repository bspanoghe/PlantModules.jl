using ModelingToolkit, OrdinaryDiffEq

t = ModelingToolkit.t_nounits
D = ModelingToolkit.D

n = 1000
@variables y(t)[1:n] = zeros(n)
@parameters jitter[1:n] = rand(n)
@variables doubledjitter(t)[1:n]


jitter_eqs = [doubledjitter[i] ~ 2*jitter[i] for i in eachindex(jitter)]
sin_eqs = [D(y[i]) ~ sin(t + doubledjitter[i]) for i in 1:length(y)]
eqs = vcat(jitter_eqs, sin_eqs)

@named sys = ODESystem(eqs, t)
sys_simpl = structural_simplify(sys);
prob = ODEProblem(sys_simpl, [], (0.0, 100))
sol = solve(prob);

@time sol[y[rand(1:n)]]
@time sol[doubledjitter[rand(1:n)]]