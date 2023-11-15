using Pkg; Pkg.activate(".")

using ModelingToolkit, DifferentialEquations, Plots, Unitful

@variables t, [description = "Time", unit = u"h"];
D = Differential(t);

# globals

function plant_compartment(; name)
    @parameters (
        R = 8.314, [description = "Ideal gas constant", unit = u"J / K * mol"], #! better way to add constant with unit?
        T = 293.15, [description = "Temperature", unit = u"K"],
        ϵ_V = 1.0, [description = "Volumetric elastic modulus", unit = u"MPa"],
        ϕ_V = 0.15, [description = "Volumetric extensibility", unit = u"MPa^-1 * hr^-1"],
        Γ = 0.3, [description = "Critical turgor pressure", unit = u"MPa"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t), [description = "Hydrostatic potential", unit = u"MPa"],
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / L"],
        ΔM(t), [description = "Change in osmotically active metabolite content", unit = u"mol / L / hr"],
        W(t), [description = "Water content", unit = u"g"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"]
    )

    eqs = [
        Ψ ~ Π + P,
        Π ~ -R*T*M,
        ΔP ~ ϵ_V / W * ΔW - ϵ_V * ϕ_V * max(P - Γ, 0),
        D(P) ~ ΔP,
        D(W) ~ ΔW,
    ]
    return ODESystem(eqs, t; name)
end

# make organs
@named root = plant_compartment()
@named stem = plant_compartment()
@named leaf = plant_compartment()

# connect via water potential
@parameters ( #! these are mostly guesses
    Ψ_soil = -0.08,
    Ψ_air = -100,
    K_soil_root = 1000,
    K_root_stem = 800,
    K_stem_leaf = 600,
    K_leaf_air = 10e-5 
)

# makes double calculations
# connections = [
#     root.F ~ K_soil_root * (Ψ_soil - root.Ψ) + K_root_stem * (stem.Ψ - root.Ψ),
#     stem.F ~ K_root_stem * (root.Ψ - stem.Ψ) + K_stem_leaf * (stem.Ψ - stem.Ψ),
#     leaf.F ~ K_soil_root * (Ψ_soil - root.Ψ) + K_root_stem * (stem.Ψ - root.Ψ),
# ]

connections = [
    F_soil2root ~ K_soil_root * (Ψ_soil - root.Ψ),
    F_root2stem ~ K_root_stem * (root.Ψ - stem.Ψ),
    F_stem2leaf ~ K_stem_leaf * (stem.Ψ - leaf.Ψ),
    F_leaf2air ~ K_leaf_air * (leaf.Ψ - Ψ_air),
    
    root.ΔW ~ F_soil2root - F_root2stem,
    stem.ΔW ~ F_root2stem - F_stem2leaf,
    leaf.ΔW ~ F_stem2leaf - F_leaf2air,

    root.ΔM ~ 1,
    stem.ΔM ~ 1,
    leaf.ΔM ~ 1
]



plant = compose(ODESystem(connections, name = :plant), root, branch, shoot1, shoot2)
plant_simp = structural_simplify(plant)

full_equations(plant_simp)

u0 = [root.X => 0.1,
      root.V => 0.2,
      branch.X => 0.1,
      branch.V => 0.2,
      shoot1.X => 0.1,
      shoot1.V => 0.2,
      shoot2.X => 0.1,
      shoot2.V => 0.2,
      ]
      
prob = ODEProblem(plant_simp, u0, (0.0, 5.0))

plot(solve(prob))