using Infiltrator
using Revise
using Plots
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful
using DataInterpolations

# # Preparatory calculations

begin
    # ## dimensions of compartments

    # ### stem segments
    segment_length = 0.2 * 1e2 # cm; from m (input)
    tree_length = 13.2 * 1e2 # cm; from m (input)
    n_segments = tree_length ÷ segment_length

    # ### stem radius at tree base
    stem_r0 = 0.2 * 1e2 / 2 # cm; from the diameter in m (input)
    radial_fraction_sapwood = 0.30 # (estimate)
    sapwood_r0 = stem_r0 * radial_fraction_sapwood

    # ### crown
    crown_base_height = 6.0 * 1e2 # cm; from m (input)
    crown_length = tree_length - crown_base_height
    n_crown_segments = crown_length ÷ segment_length
    crown_start_age = (n_segments - n_crown_segments)

    # ### branching
    branching_frequency = 3 # (estimate)
    n_branches = n_crown_segments ÷ branching_frequency

    Δr = sapwood_r0/(n_branches + 1) # (+1 to prevent stem radius from hitting 0)

    # ### branch radius
    # calculate using da vinci's rule: total cross area remains constant before and after branching
    function branch_radius(r, Δr)
        area_0 = r^2*pi
        area_new = (r - Δr)^2*pi
        area_branch = area_0 - area_new
        return sqrt(area_branch / pi)
    end

    # ### needle areas
    total_needle_area = 50.8 * 1e2^2 # cm^2; from m^2 (input)
    needle_areas = [i * total_needle_area / sum(1:n_branches) for i in n_branches:-1:1]

    # ### branchgrowth frequency
    branchgrowth_frequency = 3 # how many steps it takes branches to grow a new segment (estimate)

    # ### root dimensions
    total_root_area = 96.1 * 1e2^2 # cm^2; from m^2 (input)

    # ## parameters

    # ### hydraulic conductivities
    permeability = 2.5e-12 * 1e2^2 # cm^2; from m^2 (input)
    η = 1.0 * 1e-9 * (1/3600) # MPa hr; from mPa s (dynamic viscosity of water)
    l = segment_length # cm
    ρ_w = 1.0 # g / cm^3 (density of water)
    K_s_stem = ρ_w * permeability / η / l # g / hr / MPa / cm^2

    L_p = 3.2e-8 * 1e2 * 3600 # cm / hr / MPa; from m / s / MPa (input)
    K_roots = ρ_w * L_p * total_root_area # g / hr / MPa

    # ### elastic modulus
    ϵ_D_stem = [0.1 * 1e3, 17.5 * 0.1 * 1e3] # MPa; from GPa (input)
end

struct HollowCylinder{T} <: PlantModules.Shape
    frac_sapwood::T # fraction of radius that is conducting sapwood
    function HollowCylinder(frac_sapwood::T) where {T <: Number}
        @assert 0 <= frac_sapwood <= 1
        return new{T}(frac_sapwood)
    end
end

PlantModules.getdimensionality(::HollowCylinder) = 2
PlantModules.cross_area(hc::HollowCylinder, D::AbstractArray) = D[1]^2 * pi * (2 / hc.frac_sapwood - 1)
PlantModules.volume(hc::HollowCylinder, D::AbstractArray) = cross_area(hc, D) * D[2]

# # Structure definition

# ## Structural module types

mutable struct StemTip <: Node
    D::Vector
    height
    branch_nr
    age
end

mutable struct Stem <: Node
    D::Vector
    height
end

mutable struct BranchTip <: Node
    D::Vector
    height
    age
    needle_area
end

mutable struct Branch <: Node
    D::Vector
    height
end

# ## Growing rules

branching_rule = Rule( # top stem segment grows another thinner stem segment and a branch every `branching_frequency` steps, if stem is old enough to form crown
    StemTip,
    lhs = st -> data(st).age >= crown_start_age && data(st).age % branching_frequency == 0,
    rhs = st -> Stem(data(st).D, data(st).height) + (
        BranchTip([branch_radius(data(st).D[1], Δr), segment_length], data(st).height, 0, needle_areas[data(st).branch_nr]),
        StemTip(data(st).D - [Δr, 0.0], data(st).height + segment_length, data(st).branch_nr + 1, data(st).age + 1)
    )
)

vertical_growth_rule = Rule( # top stem segment grows only another stem segment if not branching
    StemTip,
    lhs = st -> !(data(st).age >= crown_start_age && data(st).age % branching_frequency == 0),
    rhs = st -> Stem(data(st).D, data(st).height) + 
        StemTip(data(st).D, data(st).height + segment_length, data(st).branch_nr, data(st).age + 1)
)

branchgrowth_rule = Rule( # top branch segments grown another branch segment every `branchgrowth_frequency` steps
    BranchTip,
    lhs = br -> data(br).age % branchgrowth_frequency == 0,
    rhs = br -> Branch(data(br).D, data(br).height) +
        BranchTip(data(br).D, data(br).height, data(br).age + 1, data(br).needle_area)
)

branchage_rule = Rule( # branch tips also age if not growing
    BranchTip,
    lhs = br -> !(data(br).age % branchgrowth_frequency == 0),
    rhs = br -> BranchTip(data(br).D, data(br).height, data(br).age + 1, data(br).needle_area)
)

# ## Instantiation

axiom = StemTip([sapwood_r0, segment_length], segment_length/2, 1, 0)
plant = Graph(axiom = axiom, rules = 
    (vertical_growth_rule, branching_rule, branchage_rule, branchgrowth_rule)
)
for _ in 1:(n_segments)
    rewrite!(plant)
end

PlantGraphs.prune!(plant.graph, getid([node for node in getnodes(plant) if getstructmod(node) == :StemTip][1]))
plotstructure(plant)


struct Soil <: Node end
struct Air <: Node end

graphs = [plant, Soil(), Air()]
intergraph_connections = [(1, 2) => (getnodes(plant)[1], :Soil), (1, 3) => (:BranchTip, :Air)]
plantstructure = PlantStructure(graphs, intergraph_connections)

plotstructure(plantstructure, names = "")

# # Function

# ## New functional modules

import PlantModules: t, d
readcsv(filepath) = readlines(filepath) .|> x -> split(x, ", ") .|> x -> parse(Float64, x)

# ### needle area

function needle_area_module(; name, needle_area)
    @parameters begin
        needle_area_val = needle_area, [description = "Value of total needle area on branch", unit = u"cm^2"]
    end 
    @variables begin
        needle_area(t), [description = "Total needle area on branch", unit = u"cm^2"]
    end
    eqs = [needle_area ~ needle_area_val]
    return System(eqs, t; name)
end

# ### hydraulic_module with gravity effects

function gravitational_hydraulic_module(; name, shape, ϕ_D, ϵ_D, Γ, T, D, P_0, M, height)
    D, ϕ_D, ϵ_D = [PlantModules.correctdimensionality(shape, var) for var in [D, ϕ_D, ϵ_D]] 
        # turns scalar values into vectors of correct length

    num_D = getdimensionality(shape)
    R = 8.314
    ρ_w = 1.0
    g = 9.8 * 1e-5
    P = P_0 + R*T*M - ρ_w * g * height

    @constants (
        R = R, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
        P_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / cm^3"],
        g = g, [description = "Gravitational acceleration", unit = u"hN / g"] # (from N / kg) Pa = N/m^2 => MPa = hN/cm^2
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϕ_D[1:num_D] = ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        ϵ_D[1:num_D] = ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
        height = height, [description = "Height above ground level", unit = u"cm"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        Pₕ(t), [description = "Gravitational water potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation ()
        W(t) = volume(shape, D) / ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr", guess = 0.0], #!
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P + Π + Pₕ,
        Π ~ -R*T*M,
        Pₕ ~ ρ_w * g * height,
        ΔW ~ ΣF,
        V ~ W / ρ_w,
        V ~ volume(shape, D),
        [ΔD[i] ~ D[i]*ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100) + D[i]*ΔP/ϵ_D[i] for i in eachindex(D)]...,

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return System(eqs, t; name)
end

# ### Fixed transpiration

# #### transpiration rate from data
transpiration_data = readcsv("./tutorials/clouddata/transpiration_data.csv")
transpiration_rate = LinearInterpolation(last.(transpiration_data), first.(transpiration_data))
plot(transpiration_rate) # in mg / s / m^2

# #### module
function fixed_transpiration_connection(; name, original_order)
    @constants begin
        t_unit = 1, [description = "Dummy constant for correcting units", unit = u"hr"]
        uc = 1e-3 * 3600 * (1e-2)^2, [description = "Unit conversion from (mg / s / m^2) to (g / hr / cm^2)", unit = u"g / hr / cm^2"]
    end
    @variables begin
        F_s(t), [description = "Specific water flux from compartment 2 to compartment 1", unit = u"g / hr / cm^2"]
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"]
        needle_area(t), [description = "Total needle area on branch", unit = u"cm^2"]
    end

    polarity = original_order ? -1 : 1

    eqs = [
        F_s ~ uc * transpiration_rate(t / t_unit),
        F ~ polarity * F_s * needle_area,
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, original_order) = (
        original_order ? [connection_MTK.needle_area ~ node_MTK.needle_area] : [connection_MTK.needle_area ~ nb_node_MTK.needle_area]
    )

    return System(eqs, t; name), get_connection_eqset
end


# ## Coupling

module_coupling = Dict(
	:Stem => [gravitational_hydraulic_module, constant_carbon_module, K_module],
    :Branch => [gravitational_hydraulic_module, constant_carbon_module, K_module],
	:BranchTip => [gravitational_hydraulic_module, needle_area_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module],
)

connecting_modules = Dict(
	(:Soil, :Stem) => constant_hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Branch) => hydraulic_connection,
    (:Branch, :Branch) => hydraulic_connection,
	(:Branch, :BranchTip) => hydraulic_connection,
	(:BranchTip, :Air) => fixed_transpiration_connection,
)

plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# # Parameters
begin
    module_defaults = Dict(
        :Stem => Dict(:ϵ_D => ϵ_D_stem, :K_s => K_s_stem),
        :Branch => Dict(:ϵ_D => ϵ_D_stem, :K_s => K_s_stem),
        :BranchTip => Dict(:ϵ_D => ϵ_D_stem, :K_s => K_s_stem),
        :Soil => Dict(:W_max => 1e6, :W_r => 0.5),
    )

    default_changes = Dict(
        :needle_area => 0.0, :height => 0.0, :ϕ_D => 0.0, :M => 0.0, :P_0 => PlantModules.soilfunc(0.5),
        :shape => HollowCylinder(radial_fraction_sapwood))

    connection_values = Dict(
        (:Soil, :Stem) => Dict(:K => K_roots),
    )

    plantparams = PlantParameters(; default_changes, module_defaults, connection_values)
end

# # Run it

system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false)
prob = ODEProblem(system, [], (0.0, 24.0), sparse = true)

# @time remake_graphsystem!(prob, system, plantstructure, :ϕ_D, :Stem, [10.0, 10.0])
@time sol = solve(prob, FBDF())

# # Show it

plotgraph(sol, plantstructure, structmod = [:Soil, :Stem, :Branch], varname = :Ψ)

plotgraph(sol, plantstructure, structmod = [:Stem], varname = :P, xticks = 0.0:3.0:24.0)

begin
    dimension_variables = get_subsystem_variables(system, plantstructure, :D, :Stem)

    diameter_change_mm(var, sol) = 10*(var - sol[var][1])
    low_segment_nr = ceil(Int64, 2.5 / 0.2)

    low_diameter_data = readcsv("./tutorials/clouddata/diameter_data_d.csv")
    begin
        p_diameter = plot(ylabel = "Diameter change (mm)", size = (800, 400), margins = 5*Plots.mm)#, ylims = (-0.07, 0.01), yticks = -0.07:0.01:0.01)
        plot!(p_diameter, sol, idxs = [diameter_change_mm(dimension_variables[low_segment_nr][1], sol)], label = "Simulated")
        plot!(p_diameter, first.(low_diameter_data), last.(low_diameter_data), label = "Data")
    end
end

plot(sol, idxs = [get_subsystem_variables(system, plantstructure, :F, (:Stem, :Stem))[2] / 1000],
    xlims = (6.0, 18.0), ylabel = "Sap flow (kg/h)", xlabel = "Time of day (h)", legend = false)

plotgraph(sol, plantstructure, structmod = :Soil, varname = :ΣF)
plotgraph(sol, plantstructure, structmod = :Soil, varname = :Ψ)

plotgraph(sol, plantstructure, structmod = :Air, varname = :ΣF)

# check transpiration
transp_max_paper = 8.0 # mg / m^2 / s
transp_max_expected = transp_max_paper * 1e-3 * (total_needle_area * (1e-2)^2) * 3600 # g / hr