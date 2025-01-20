using Pkg; Pkg.activate("./tutorials")
using PlantModules
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as d

@mtkmodel mini_hydraulics begin
    @parameters begin
        T = 298.15, [description = "Temperature"]
        ϕ_D = 3.0, [description = "extensibility"]
        Γ = 0.3, [description = "critical pressure"]
    end
    @variables begin
        P(t), [description = "Hydrostatic potential"]
        W(t), [description = "Water content"]
        D(t), [description = "Dimension of compartment"]
        ΣF(t), [description = "Net incoming water flux"]
        
        ΔP(t), [description = "Change in hydrostatic potential"]
        ΔW(t), [description = "Change in water content"]
        ΔD(t), [description = "Change in dimensions of compartment"]
    end

    @equations begin
        ΔW ~ ΣF
        W ~ D^3
        ΔD ~ D * (ΔP + ϕ_D * (P - Γ))
        
        d(P) ~ ΔP
        d(W) ~ ΔW
        d(D) ~ ΔD
    end
end

@mtkmodel mini_connection begin
    @parameters begin
        K = K, [description = "Hydraulic conductivity of connection"]
    end
    @variables begin
        F(t), [description = "Water flux from compartment 1 to compartment 2"]
        P_1(t), [description = "Total water potential of compartment 1"]
        P_2(t), [description = "Total water potential of compartment 2"]
    end
    @equations begin
        F ~ K * (P_1 - P_2)
    end
end

# Create compartment instances #

@named stem = mini_hydraulics(P = 0.1, D = 3.0)
@named leaf = mini_hydraulics(P = 0.1, D = 1.0)

@named stem_leaf = mini_connection(K = 600)
# @named leaf_stem = mini_connection(K = 600)

# define connections #

connections = [
    stem_leaf.P_1 ~ stem.P,
    stem_leaf.P_2 ~ leaf.P,
    stem.ΣF ~ - stem_leaf.F,

    # leaf_stem.P_1 ~ leaf.P,
    # leaf_stem.P_2 ~ stem.P,
    # leaf.ΣF ~ - leaf_stem.F,
    leaf.ΣF ~ stem_leaf.F
]

# build model #

## model definition
plant = compose(ODESystem(connections, t, name = :plant), stem, leaf, stem_leaf)#, leaf_stem)
plant_simp = structural_simplify(plant)

# full_equations(plant_simp)

## define and solve problem
prob = ODEProblem(plant_simp, [], (0.0, 24.0*31), guesses = [stem.ΔP => 0.0, leaf.ΔP => 0.0])
sol = solve(prob)

plot(sol)