using Pkg; Pkg.activate(".")

using ModelingToolkit, DifferentialEquations, Plots, Unitful

@variables t, [description = "Time", unit = u"hr"];
d = Differential(t);

# define plant compartment volumes #

abstract type Volume end

# Write dimensions in the order: radius
struct Sphere<:Volume
    ϵ_D::Vector
    ϕ_D::Vector
    formula::Function
    function Sphere(ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 1 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 1 was expected.")
        length(ϕ_D) != 1 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 1 was expected.")

        return new(ϵ_D, ϕ_D, D -> 4/3 * pi * D[1]^3)
    end
end

# Write dimensions in the order: radius - length
struct Cilinder<:Volume
    ϵ_D::Vector
    ϕ_D::Vector
    formula::Function
    function Cilinder(ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 2 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 2 was expected.")
        length(ϕ_D) != 2 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 2 was expected.")

        return new(ϵ_D, ϕ_D, D -> D[1]^2 * pi * D[2])
    end
end

# Write dimensions in the order: length - width - height
struct Rectangle<:Volume
    ϵ_D::Vector
    ϕ_D::Vector
    formula::Function
    function Rectangle(ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 3 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 3 was expected.")
        length(ϕ_D) != 3 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 3 was expected.")

        return new(ϵ_D, ϕ_D, D -> D[1] * D[2] * D[3])
    end
end

# define plant compartments #

## constants
@constants (
    R = 8.314e-6, [description = "Ideal gas constant", unit = u"MJ / K / mol"],
    ρ_w = 1.0e6,  [description = "Density of water", unit = u"g / m^3"]
)

## general compartment definition
function plant_compartment(; name, volume::Volume)
    num_D = length(volume.ϵ_D)
    @parameters (
        T = 298.15, [description = "Temperature", unit = u"K"],
        # ϵ_V = 1.0, [description = "Volumetric elastic modulus", unit = u"MPa"],
        ϵ_D[1:num_D] = volume.ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        # ϕ_V = 0.15, [description = "Volumetric extensibility", unit = u"MPa^-1 * hr^-1"],
        ϕ_D[1:num_D] = volume.ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        Γ = 0.3, [description = "Critical turgor pressure", unit = u"MPa"],
        P_0 = 0.0, [description = "Minimum pressure", unit = u"MPa"], 
            # you can't use magic numbers in the equations when working with units
            # everything has to be either a variable or a parameter
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t), [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / m^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
        W(t), [description = "Water content", unit = u"g"],
        D(t)[1:num_D], [description = "Dimensions of compartment", unit = u"m"],
        V(t), [description = "Volume of compartment", unit = u"m^3"],

        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        ΔM(t), [description = "Change in osmotically active metabolite content", unit = u"mol / m^3 / hr"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"m / hr"],
        ΔV(t), [description = "Change in volume of compartment", unit = u"m^3 / hr"],

        F(t), [description = "Net incoming water flux", unit = u"g / hr"],
    )

    eqs = [
        Ψ ~ Π + P, # Water potential consists of a solute- and a pressure component
        Π ~ -R*T*M, # Solute component is determined by concentration of dissolved metabolites
        [ΔP ~ ϵ_D[i] * (ΔD[i]/D[i] - ϕ_D[i]*max(P - Γ, P_0)) for i in eachindex(D)]..., # Pressure changes with changes in compartment dimensions
            #! im assuming the max is also required here
        V ~ volume.formula(D), # Compartment dimensions determine its volume
        W ~ ρ_w * V, # Volume determines water content (*ignoring dry matter)
        ΔW ~ F, # Water content changes due to water fluxes, determined by water potentials

        d(P) ~ ΔP,
        d(M) ~ ΔM,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
        d(W) ~ ΔW,
        d(V) ~ ΔV
    ]
    return ODESystem(eqs, t; name)
end

## named compartments
@named root = plant_compartment(volume = Sphere([3.0], [0.45]))
@named stem = plant_compartment(volume = Cilinder([6.0, 0.15], [0.8, 0.03]))
@named leaf = plant_compartment(volume = Rectangle([5.0, 0.3, 0.2], [0.7, 0.1, 0.05]))


# define connections #


## connection parameters
### constant parameters
@parameters ( #! these are mostly guesses
    Ψ_soil = -0.3, [description = "Total water potential of soil", unit = u"MPa"],
    Ψ_air = -100, [description = "Total water potential of air", unit = u"MPa"],
    K_soil_root = 50, [description = "Hydraulic conductivity between soil and root", unit = u"g / hr / MPa"],
    K_root_stem = 800, [description = "Hydraulic conductivity between root and stem", unit = u"g / hr / MPa"],
    K_stem_leaf = 600, [description = "Hydraulic conductivity between stem and leaf", unit = u"g / hr / MPa"],
    K_leaf_air = 1e-3, [description = "Hydraulic conductivity between leaf and air", unit = u"g / hr / MPa"],
    A_max = 10, [description = "Maximum rate of photosynthesis", unit = u"mol / m^3 / hr"],
    A_0 = 0, [description = "Rate of photosynthesis if there is no photosynthesis", unit = u"mol / m^3 / hr"],
    R = 1.5, [description = "Rate of cellular respiration", unit = u"mol / m^3 / hr"],
    K_M = 0.1, [description = "Rate of metabolite diffusion", unit = u"hr^-1"],
)

### variable parameters
val(x) = x
val(x::Quantity) = x.val

A_n(t, A_max) = A_max/2 * (sin(val(t) * pi/12 - pi/2) + 1) # simulate day and night cycle of light
# plot(t -> A_n(t, 5), xlims = (0, 48), xticks = [0, 12, 24, 36, 48])

@register_symbolic A_n(t, A_max)

## connections themselves
connections = [
    root.F ~ K_soil_root * (Ψ_soil - root.Ψ) + K_root_stem * (stem.Ψ - root.Ψ),
    stem.F ~ K_root_stem * (root.Ψ - stem.Ψ) + K_stem_leaf * (leaf.Ψ - stem.Ψ),
    leaf.F ~ K_stem_leaf * (stem.Ψ - leaf.Ψ) + K_leaf_air * (Ψ_air - leaf.Ψ),

    root.ΔM ~ A_0 - R + K_M * (stem.M - root.M),
    stem.ΔM ~ A_0 - R + K_M * (root.M - stem.M) + K_M * (leaf.M - stem.M),
    leaf.ΔM ~ A_n(t, A_max) - R + K_M * (stem.M - leaf.M)
]


# build model #

## model definition
plant = compose(ODESystem(connections, name = :plant), root, stem, leaf)
plant_simp = structural_simplify(plant)

# full_equations(plant_simp)

## initial values
u0 = [
    root.M => 200.0,
    root.P => 0.1,
    root.W => 3.0,
    root.D[1] => 0.1,

    stem.M => 200.0,
    stem.P => 0.1,
    stem.W => 3.0,
    stem.D[1] => 0.4,
    stem.D[2] => 0.03,

    leaf.M => 200.0,
    leaf.P => 0.1,
    leaf.W => 3.0,
    leaf.D[1] => 0.3,
    leaf.D[2] => 0.05,
    leaf.D[3] => 0.03
]

### Adding initial values for dummy derivatives generated by MTK (see https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/#ERROR:-ArgumentError:-SymbolicUtils.BasicSymbolic{Real}[x%CB%8Dt(t)]-are-missing-from-the-variable-map.)
u0_req = ModelingToolkit.missing_variable_defaults(plant_simp) # generates zero initial value for ALL initial states
u0_ext = union(u0, u0_req) # add actual values of non-dummy initial states 
u0_full = unique(x -> x[1], u0_ext) # remove duplicates

## define and solve problem
prob = ODEProblem(plant_simp, u0_full, (0.0, 24.0*31))
sol = solve(prob)

## plot solution
plots = [plot(sol, idxs = [getproperty(organ, var)]) for organ in [root, stem, leaf] for var in [:W, :P, :M, :Ψ]]
plot(plots..., layout = (3, 4), size = (800, 500))

D_plots = [plot(sol, idxs = [organ.D...]) for organ in [root, stem, leaf]]
plot(D_plots..., layout = (3, 1), size = (800, 500))

#=
# other approach to defining connections #

function plant_connection(; name, c1, c2, K)
    @parameters (
        K = K, [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 1 to compartment 2", unit = u"g / hr"],
    )

    eqs = [
        F ~ K * (c1.Ψ - c2.Ψ)
    ]
    return ODESystem(eqs, t; name)
end

@named soil_root = plant_connection(c1 = root, c2 = stem, K = 1000)
@named root_stem = plant_connection(K = 800)
@named stem_leaf = plant_connection(K = 600)
@named leaf_air = plant_connection(K = 1e-5)

connections = [
    soil_root.F ~ K_soil_root * (Ψ_soil - root.Ψ),
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
=#