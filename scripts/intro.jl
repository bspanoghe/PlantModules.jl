using Pkg; Pkg.activate(".")

using ModelingToolkit, DifferentialEquations, Plots, Unitful

@variables t, [description = "Time", unit = u"hr"];
d = Differential(t);

# Define plant compartment shapes #

abstract type Shape end

struct Sphere<:Shape
    ϵ_D::Vector
    ϕ_D::Vector
    function Sphere(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 1 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 1 was expected.")
        length(ϕ_D) != 1 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 1 was expected.")

        return new(ϵ_D, ϕ_D)
    end
end

struct Cilinder<:Shape
    ϵ_D::Vector
    ϕ_D::Vector
    function Cilinder(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 2 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 2 was expected.")
        length(ϕ_D) != 2 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 2 was expected.")

        return new(ϵ_D, ϕ_D)
    end
end

struct Cuboid<:Shape
    ϵ_D::Vector
    ϕ_D::Vector
    function Cuboid(; ϵ_D::Vector, ϕ_D::Vector)
        length(ϵ_D) != 3 && error("An array of length $(length(ϵ_D)) was given for ϵ_D while length 3 was expected.")
        length(ϕ_D) != 3 && error("An array of length $(length(ϕ_D)) was given for ϕ_D while length 3 was expected.")

        return new(ϵ_D, ϕ_D)
    end
end

volume(s::Shape, D::AbstractArray) = error("Function volume is not defined for shape $s")
volume(::Sphere, D::AbstractArray) = 4/3 * pi * D[1]^3 # Write dimensions in the order: radius
volume(::Cilinder, D::AbstractArray) = D[1]^2 * pi * D[2] # Write dimensions in the order: radius - length
volume(::Cuboid, D::AbstractArray) = D[1] * D[2] * D[3] # Write dimensions in the order: length - width - height

# Define helper functions #

## Unitful is a dangerous beast
val(x) = x
val(x::Quantity) = x.val

## Forbidden rites
import Base.exp
import Base.log
exp(x::Quantity) = exp(val(x))*unit(x)
log(x::Quantity) = log(val(x))*unit(x)

"""
LSE(x, y, ...; γ = 1)

LogSumExp, a smooth approximation for the maximum function. 
The temperature parameter γ determines how close the approximation is to the actual maximum.
WARNING: large values for γ combined with large input arguments
will result in numerical overflow and cause the function to return Inf.
"""
LSE(x::Real...; γ = 1) = log(sum(exp.(γ .* x)) ) / γ 
    #! add warning for users in case pressures are incorrect unit?
@register_symbolic LSE(x)

# Define compartments #

## Plants

### constants

@constants (
    R = 8.314e-6, [description = "Ideal gas constant", unit = u"MJ / K / mol"],
)

### general compartment definition
function plant_compartment(; name, shape::Shape)
    num_D = length(shape.ϵ_D)
    @parameters (
        T = 298.15, [description = "Temperature", unit = u"K"],
        ρ_w = 1.0e6, [description = "Density of water", unit = u"g / m^3"],
        ϵ_D[1:num_D] = shape.ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        ϕ_D[1:num_D] = shape.ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        Γ = 0.3, [description = "Critical turgor pressure", unit = u"MPa"],
        P_0 = 0.0, [description = "Minimum pressure", unit = u"MPa"], 
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t), [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / m^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
        W(t), [description = "Water content", unit = u"g"],
        D(t)[1:num_D], [description = "Dimensions of compartment", unit = u"m"],
        V(t), [description = "Shape of compartment", unit = u"m^3"],
        F(t), [description = "Net incoming water flux", unit = u"g / hr"],

        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        ΔM(t), [description = "Change in osmotically active metabolite content", unit = u"mol / m^3 / hr"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"m / hr"],
    )

    eqs = [
        Ψ ~ Π + P, # Water potential consists of a solute- and a pressure component
        Π ~ -R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ F, # Water content changes due to flux (depending on water potentials as defined in connections)
        V ~ W / ρ_w, # Shape is directly related to water content        
        V ~ volume(shape, D), # Shape is also directly related to compartment dimensions
        [ΔD[i] ~ D[i] * (ΔP/ϵ_D[i] + ϕ_D[i] * LSE(P - Γ, P_0, γ = 100)) for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure

        d(P) ~ ΔP,
        d(M) ~ ΔM, # Change in dissolved metabolites is defined in the connections
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return ODESystem(eqs, t; name)
end

## Environmental compartments

function environmental_compartment(; name, W_max)
    @parameters (
        T = 298.15, [description = "Temperature", unit = u"K"],
        ρ_w = 1.0e6, [description = "Density of water", unit = u"g / m^3"],
        W_max = W_max, [description = "Water capacity of compartment", unit = u"g"],
        )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        W(t), [description = "Water content", unit = u"g"],
        W_r(t), [description = "Relative water content", unit = u"g / g"],
        F(t), [description = "Net incoming water flux", unit = u"g / hr"],

        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
    )

    eqs = [
        W_r ~ W / W_max,
        ΔW ~ F, # Water content changes due to flux (depending on water potentials as defined in connections)
        d(W) ~ ΔW,

        # prevent variables from being erased from existence
        Ψ ~ Ψ,
        T ~ T
    ]
    return ODESystem(eqs, t; name)
end

## Compartment connections

function compartment_connection(; name, K)
    @parameters (
        K = K, [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 1 to compartment 2", unit = u"g / hr"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ K * (Ψ_1 - Ψ_2)
    ]
    return ODESystem(eqs, t; name)
end

# Create compartment instances #

rootvol = Sphere(ϵ_D = [3.0], ϕ_D = [0.45])
stemvol = Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.8, 0.03])
leafvol = Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.7, 0.1, 0.05])

@named root = plant_compartment(shape = rootvol)
@named stem = plant_compartment(shape = stemvol)
@named leaf = plant_compartment(shape = leafvol)

@named soil = environmental_compartment(W_max = 500) # 500 g of water max - we're growing this plant in a pot!
@named air = environmental_compartment(W_max = 1e4) # how much water can a (hopefully ventilated) house full of air hold?...

@named soil_root = compartment_connection(K = 50)
@named root_stem = compartment_connection(K = 800)
@named stem_leaf = compartment_connection(K = 600)
@named leaf_air = compartment_connection(K = 1e-3)

# define connections #


## connection parameters
### constant parameters

@parameters ( #! these are mostly guesses
    V_w = 18e-6, [description = "Molar volume of water", unit = u"m^3 / mol"],
    A_max = 10, [description = "Maximum rate of photosynthesis", unit = u"mol / m^3 / hr"],
    A_0 = 0, [description = "Rate of photosynthesis if there is no photosynthesis", unit = u"mol / m^3 / hr"],
    Rsp = 1.5, [description = "Rate of cellular respiration", unit = u"mol / m^3 / hr"],
    K_M = 0.3, [description = "Rate of metabolite diffusion", unit = u"hr^-1"],

    unit_MPa = 1.0, [description = "A dummy parameter to correct units in an equation", unit = u"MPa"],
)

### variable parameters

A_n(t, A_max) = A_max/2 * (sin(val(t) * pi/12 - pi/2) + 1) # simulate day and night cycle of light
# plot(t -> A_n(t, 5), xlims = (0, 48), xticks = [0, 12, 24, 36, 48])

@register_symbolic A_n(t, A_max)

Ψ_soil_func(W_r) = -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19) # based on https://www.researchgate.net/figure/Relationship-between-soil-water-potential-and-soil-water-content-of-field-capacity_fig1_8888669
    # adapted with 1/x to make Ψ → -Inf as W_r → 0
# plot(Ψ_soil_func, xlims = (0, 1), ylims = (-10, 0))

@register_symbolic Ψ_soil_func(W_r)

## connections themselves

connections = [
    soil.F ~ 0 - soil_root.F,
    root.F ~ soil_root.F - root_stem.F,
    stem.F ~ root_stem.F - stem_leaf.F,
    leaf.F ~ stem_leaf.F - leaf_air.F,
    air.F ~ leaf_air.F - 0,

    soil_root.Ψ_1 ~ soil.Ψ,
    soil_root.Ψ_2 ~ root.Ψ,
    root_stem.Ψ_1 ~ root.Ψ,
    root_stem.Ψ_2 ~ stem.Ψ,
    stem_leaf.Ψ_1 ~ stem.Ψ,
    stem_leaf.Ψ_2 ~ leaf.Ψ,
    leaf_air.Ψ_1 ~ leaf.Ψ,
    leaf_air.Ψ_2 ~ air.Ψ,

    root.ΔM ~ A_0 - Rsp + K_M * (stem.M - root.M),
    stem.ΔM ~ A_0 - Rsp + K_M * (root.M - stem.M) + K_M * (leaf.M - stem.M),
    leaf.ΔM ~ A_n(t, A_max) - Rsp + K_M * (stem.M - leaf.M),

    soil.Ψ ~ Ψ_soil_func(soil.W_r) * unit_MPa, #! add (a) universal soil function as default? e.g. van genuchten vergelijking
        #! is there some way to have unitful ignore an equation?
    air.Ψ ~ R * air.T / V_w * log(air.W_r) #! add multiple blocks for soil and air with infinite W_max as default
]

# build model #

## model definition
plant = compose(ODESystem(connections, name = :plant), soil, root, stem, leaf, air, soil_root, root_stem, stem_leaf, leaf_air)
plant_simp = structural_simplify(plant)

# full_equations(plant_simp)

## initial values

u0 = [
    soil.W => soil.W_max/2,

    root.M => 200.0,
    root.P => 0.1,
    root.D[1] => 0.1,
    root.W => volume(rootvol, root.D) * root.ρ_w,

    stem.M => 200.0,
    stem.P => 0.1,
    stem.D[1] => 0.4,
    stem.D[2] => 0.03,
    stem.W => volume(stemvol, stem.D) * stem.ρ_w,

    leaf.M => 200.0,
    leaf.P => 0.1,
    leaf.D[1] => 0.3,
    leaf.D[2] => 0.05,
    leaf.D[3] => 0.03,
    leaf.W => volume(leafvol, leaf.D) * leaf.ρ_w,

    air.W => air.W_max/2,
]

### Adding initial values for dummy derivatives generated by MTK (see https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/#ERROR:-ArgumentError:-SymbolicUtils.BasicSymbolic{Real}[x%CB%8Dt(t)]-are-missing-from-the-variable-map.)
u0_req = ModelingToolkit.missing_variable_defaults(plant_simp) # generates zero initial value for ALL initial states
u0_ext = union(u0, u0_req) # add actual values of non-dummy initial states 
u0_full = unique(x -> x[1], u0_ext) # remove duplicates

## define and solve problem
prob = ODEProblem(plant_simp, u0_full, (0.0, 24.0*31))
sol = solve(prob)

## plot solution
org_plots = [plot(sol, idxs = [getproperty(organ, var)]) for organ in [root, stem, leaf] for var in [:W, :P, :M, :Ψ]]
plot(org_plots..., layout = (3, 4), size = (800, 500))

env_plots = [plot(sol, idxs = [getproperty(organ, var)]) for organ in [soil, root, stem, leaf, air] for var in [:W, :Ψ]]
plot(env_plots..., layout = (5, 2), size = (800, 500))

D_plots = [plot(sol, idxs = [organ.D...]) for organ in [root, stem, leaf]]
plot(D_plots..., layout = (3, 1), size = (800, 500))

flux_plots = [plot(sol, idxs = [getproperty(connection, :F)]) for connection in [soil_root, root_stem, stem_leaf, leaf_air]]
plot(flux_plots..., layout = (4, 1), size = (800, 500))

### repeat for larger timescale

prob = ODEProblem(plant_simp, u0_full, (0.0, 24.0*365*20)) # i put my plant in a pot and leave it there for 20 years
sol = solve(prob)

env_plots = [plot(sol, idxs = [getproperty(organ, var)]) for organ in [soil, root, stem, leaf, air] for var in [:W, :Ψ]]
plot(env_plots..., layout = (5, 2), size = (800, 500))