using Pkg; Pkg.activate("./tutorials")
using PlantModules
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as d

function mini_hydraulics(; name, D)
    @variables (
        D(t) = D, [description = "Water height"],
        ΣF(t), [description = "Net incoming water flux"],
    )

    eqs = [
        d(D) ~ ΣF
    ]
    return ODESystem(eqs, t; name)
end

function mini_connection(; name, K)
    @parameters (
        K = K, [description = "Hydraulic conductivity of connection"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 1 to compartment 2"],
        D_1(t), [description = "Total water potential of compartment 1"],
        D_2(t), [description = "Total water potential of compartment 2"],
    )

    eqs = [
        F ~ K * (D_1 - D_2)
    ]
    return ODESystem(eqs, t; name)
end

# Create compartment instances #

@named stem = mini_hydraulics(D = 3.0)
@named leaf = mini_hydraulics(D = 0.5)

@named stem_leaf = mini_connection(K = 6)
@named leaf_stem = mini_connection(K = 6)

# define connections #

connections = [
    stem_leaf.D_1 ~ stem.D,
    stem_leaf.D_2 ~ leaf.D,
    stem.ΣF ~ - stem_leaf.F,

    leaf_stem.D_1 ~ leaf.D,
    leaf_stem.D_2 ~ stem.D,
    leaf.ΣF ~ - leaf_stem.F,
]

# build model #

## model definition
plant = compose(ODESystem(connections, t, name = :plant), stem, leaf, stem_leaf, leaf_stem)
plant_simp = structural_simplify(plant)

# full_equations(plant_simp)

## define and solve problem
prob = ODEProblem(plant_simp, [], (0.0, 10.0))
sol = solve(prob)

plot(sol)