using Revise
using Plots
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful
using DataInterpolations

# # Structure definition

# ## Preparatory calculations

# ### stem segments and change in radius
segment_length = 20 # cm input
tree_length = 1200 # cm input
stem_r0 = 12.0 # cm input
branching_frequency = 3 # dealer's choice

n_segments = tree_length ÷ segment_length
Δr = stem_r0/n_segments

# ### crown
crown_length = tree_length - 600 # cm input
n_crown_segments = crown_length ÷ segment_length
segments_before_crown = (n_segments - n_crown_segments)

# ### branch radius
function branch_radius(r, Δr)
    area_0 = r^2*pi
    area_new = (r - Δr)^2*pi
    area_branch = area_0 - area_new  # from da vinci's rule (total cross area remains constant)
    return sqrt(area_branch / pi)
end

# ### total needle length

get_cylinder_length(area, radius) = (area - 2*radius^2*pi) / (2*radius*pi)
    # A = 2 * (D[1]^2 * pi) + (2 * D[1] * pi) * D[2]
    # D[2] = (A - 2 * (D[1]^2 * pi)) /  (2 * D[1] * pi)

needle_area = 51 * 100^2 # cm^2 input
needle_radius = 0.06 # cm based on https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.13465
total_needle_length = get_cylinder_length(needle_area, needle_radius)

# ### sanity check: estimate number of needles in tree assuming 10 cm needles
num_needles = needle_area / (needle_radius^2 * pi * 10)

# ### needles
n_needle_segments = n_crown_segments ÷ 3
needle_lengths = [i * total_needle_length / sum(1:n_needle_segments) for i in n_needle_segments:-1:1]

# ### root dimensions

rootarea = 96.1 * 100^2 # cm^2 input
rootradius = 0.015 # cm based on https://cdnsciencepub.com/doi/abs/10.1139/x98-206
total_root_length = get_cylinder_length(rootarea, rootradius) #! hm.


# ## Actual definition

Base.@kwdef mutable struct StemTip <: Node
    tsb::Integer = 1 # time since branching
    nd::Integer = 1 # number of divisions / branchings
    D::Vector = [stem_r0, segment_length]
end

Base.@kwdef mutable struct Stem <: Node
    D::Vector
end

Base.@kwdef mutable struct Branch <: Node
    D::Vector
    age::Integer = 1
end

Base.@kwdef mutable struct Needles <: Node
    D::Vector
end

Base.@kwdef mutable struct Roots <: Node
    D::Vector = [rootradius, total_root_length] 
end

# ## Growing rules

# ### Plant parts

stop_rule = Rule(
    StemTip,
    lhs = gt -> data(gt).D[1] < Δr
)

vertical_growth_rule = Rule(
    StemTip,
    lhs = gt -> data(gt).tsb < (segments_before_crown+branching_frequency) && data(gt).D[1] > Δr,
    rhs = gt -> Stem(data(gt).D) + StemTip(data(gt).tsb + 1, data(gt).nd, data(gt).D - [Δr, 0.0])
)

branching_rule = Rule(
    StemTip,
    lhs = gt -> data(gt).tsb == (segments_before_crown+branching_frequency),
    rhs = gt -> Stem(data(gt).D) + (
        Branch([branch_radius(data(gt).D[1], Δr), segment_length], 0) + Needles([needle_radius, needle_lengths[data(gt).nd]]),
        StemTip(segments_before_crown+1, data(gt).nd + 1, data(gt).D - [Δr, 0.0])
    )
)

branchgrowth_rule = Rule(
    Branch,
    lhs = br -> (data(br).age == 1 && length(children(br)) == 1),
    rhs = br -> Branch(data(br).D, data(br).age + 1) + Branch(data(br).D, 0)
)

# ## Environment

struct Soil <: Node end
struct Air <: Node end

# ## Instantiation

axiom = Roots() + StemTip()
plant = Graph(axiom = axiom, rules = (stop_rule, vertical_growth_rule, branching_rule, branchgrowth_rule))
for _ in 1:(n_segments+10)
    rewrite!(plant)
end

graphs = [plant, Soil(), Air()]
intergraph_connections = [(1, 2) => (:Roots, :Soil), (1, 3) => (:Needles, :Air)]
plantstructure = PlantStructure(graphs, intergraph_connections)

plotstructure(plant)
plotstructure(plantstructure)

# # Function

# ## New functional modules

transpiration_data = readlines("./tutorials/transpiration_data.csv") .|> x -> split(x, ", ") .|> x -> parse(Float64, x)
transpiration_rate = LinearInterpolation(last.(transpiration_data), first.(transpiration_data))
plot(transpiration_rate) # in mg / s / m^2

function steered_waterpotential_module(; name, shape, ϕ_D, ϵ_D, Γ, T, D, Ψ, M)
    D, ϕ_D, ϵ_D = [PlantModules.correctdimensionality(shape, var) for var in [D, ϕ_D, ϵ_D]] 
        # turns scalar values into vectors of correct length

    num_D = getdimensionality(shape)
    R = 8.314
    P = Ψ + M*R*T

    @constants (
        R = R, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
        P_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"],
        ρ_w = 1.0, [description = "Density of water", unit = u"g / cm^3"],
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϕ_D[1:num_D] = ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        ϵ_D[1:num_D] = ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t), [description = "Hydrostatic potential", unit = u"MPa", guess = 0.0], #! initial value
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3", guess = 0.0], # m^3 so units match in second equation ()
        W(t) = volume(shape, D) / ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔM(t), [description = "Change in metabolite content", unit = u"mol / cm^3 / hr"],
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr", guess = 0.0], #!
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P + Π, # Water potential consists of a solute- and a pressure component
        Π ~ -R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        V ~ W / ρ_w, # Volume is directly related to water content  
        V ~ volume(shape, D), # Volume is also directly related to compartment dimensions
        [ΔD[i] ~ D[i]*ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100) + D[i]*ΔP/ϵ_D[i] for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure

        d(M) ~ ΔM,
        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return System(eqs, t; name)
end

function unknown_Ψ_module(; name)
    @variables begin
        Ψ(t), [description = "Total water potential", unit = u"MPa"]
    end
    eqs = [Ψ ~ Ψ]
    return System(eqs, t; name)
end

function fixed_transpiration_connection(; name, shape, K)
    num_D = getdimensionality(shape)

    @constants begin
        t_unit = 1, [description = "Dummy constant for correcting units", unit = u"hr"]
        uc = 1e-3 * 3600 * (1e-2)^2, [description = "Unit conversion from (mg / s / m^2) to (g / hr / cm^2)", unit = u"g / hr / cm^2"]
    end
    @parameters begin
        K(t) = K, [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"]
    end
    @variables begin
        F_s(t), [description = "Specific water flux from compartment 2 to compartment 1", unit = u"g / hr / cm^2"]
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"]
        D(t)[1:num_D], [description = "Dimensions of compartment", unit = u"cm"]
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"]
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"]
    end

    eqs = [
        F_s ~ uc * transpiration_rate(t / t_unit),
        F ~ F_s * surface_area(shape, D) / 2,
        F ~ K * (Ψ_2 - Ψ_1)
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, correct_order) = (
        correct_order ? 
        [
            connection_MTK.Ψ_1 ~ node_MTK.Ψ,
            connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,

            connection_MTK.D ~ node_MTK.D,
        ] : [
            connection_MTK.Ψ_1 ~ node_MTK.Ψ,
            connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,

            connection_MTK.D ~ nb_node_MTK.D,
        ]
    )

    return System(eqs, t; name), get_connection_eqset
end


# ## Coupling

module_coupling = Dict(
	:Roots => [hydraulic_module, constant_carbon_module, K_module],
	:Stem => [hydraulic_module, constant_carbon_module, K_module],
    :Branch => [hydraulic_module, constant_carbon_module, K_module],
	:Needles => [steered_waterpotential_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module],
)

connecting_modules = Dict(
	(:Soil, :Roots) => constant_hydraulic_connection,
	(:Roots, :Stem) => constant_hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Branch) => hydraulic_connection,
    (:Branch, :Branch) => hydraulic_connection,
	(:Branch, :Needles) => constant_hydraulic_connection,
	(:Needles, :Air) => fixed_transpiration_connection,
)

plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# # Parameters
begin
    default_changes = Dict{Symbol, Any}(:ϕ_D => 0.0)

    permeability = 2e-12 * 1e4 # cm^2 (from m^2)
    η = 1.0 * 1e-9 * (1/3600) # MPa hr (from mPa s)
    l = segment_length # cm
    ρ_w = 1.0 # g / cm^3
    K_s_stem = ρ_w * permeability / η / l # g / hr / MPa / cm^2

    module_defaults = Dict(
        :Stem => Dict(:ϵ_D => [0.1 * 1e3, 17.5 * 0.1 * 1e3], :K_s => K_s_stem),
        :Needles => Dict(:M => 450e-6, :area_func => surface_area, :M_c => 0.4, :A_max => 1500, :t_sunrise => 6.0, :t_sunset => 20.0),
        :Soil => Dict(:W_max => 1e8, :T => 288.15, :W_r => 0.5),
    )

    L_p = 3.2e-8 * 1e2 * 3600 # cm / hr / MPa (from m / s / MPa)
    K_roots = ρ_w * L_p * rootarea # g / hr / MPa

    connection_values = Dict(
        (:Soil, :Roots) => Dict(:K => K_roots),
        (:Roots, :Stem) => Dict(:K => 1e9), # roots => stem not modelled in original system
        (:Branch, :Needles) => Dict(:K => 1e5),
        (:Needles, :Air) => Dict(:K => 1e-3),
    )

    plantparams = PlantParameters(; default_changes, module_defaults, connection_values)
end

# # Run it

system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false)
prob = ODEProblem(system, [], (0.0, 24.0), sparse = true)
@time sol = solve(prob);

# # Show it
plotgraph(sol, plantstructure, structmod = [:Soil, :Roots, :Stem, :Branch, :Needles], varname = :Ψ)

plotgraph(sol, plantstructure, structmod = :Stem, varname = :D)
plotgraph(sol, plantstructure, structmod = :Stem, varname = :ΣF)
plotnode(sol, getnodes(plantstructure)[10], varname = :D, xlims = (10*24.0, 11*24.0), ylims = (10.39, 10.4))
plotgraph(sol, plantstructure, structmod = [:Roots, :Stem, :Branch, :Needles], varname = :W)

plotgraph(sol, plantstructure, structmod = :Needles, varname = :W)
plotgraph(sol, plantstructure, structmod = :Needles, varname = :M)
plotgraph(sol, plantstructure, structmod = :Needles, varname = :Ψ)

plotgraph(sol, plantstructure, structmod = :Soil, varname = :W)
plotgraph(sol, plantstructure, structmod = :Soil, varname = :Ψ)

plotgraph(sol, plantstructure, structmod = :Needles, varname = :ΣF)
plotgraph(sol, plantstructure, structmod = :Air, varname = :ΣF)