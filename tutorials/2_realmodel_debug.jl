using Plots
using Revise, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq
using DataInterpolations

# # Preparatory calculations

begin
    # dimensions of compartments

    ## stem segments
    segment_length = 0.2 * 1e2 # cm; from m (input)
    tree_length = 13.2 * 1e2 # cm; from m (input)
    n_segments = tree_length ÷ segment_length

    ## stem radius at tree base
    stem_r0 = 0.2 * 1e2 / 2 # cm; from the diameter in m (input)
    radial_fraction_sapwood = 0.25 # (estimate)
    sapwood_r0 = stem_r0 * radial_fraction_sapwood

    ## crown
    crown_base_height = 6.0 * 1e2 # cm; from m (input)
    crown_length = tree_length - crown_base_height
    n_crown_segments = crown_length ÷ segment_length
    crown_start_age = (n_segments - n_crown_segments)

    ## branching
    branching_frequency = 3 # (estimate)
    n_branches = n_crown_segments ÷ branching_frequency

    Δr = sapwood_r0/(n_branches + 1) 
		# how much the radius of a stem segment decreases whenever the stem branches
		# we assume the tree decreases linearly in radius from the given starting radius at the bottom to a radius of near 0 at the top
		# (+1 to prevent stem radius from hitting 0)

    ## branch radius
    	# calculate using da vinci's rule: total cross area remains constant
		# before and after branching
    function branch_radius(r, Δr)
        area_0 = r^2*pi
        area_new = (r - Δr)^2*pi
        area_branch = area_0 - area_new
        return sqrt(area_branch / pi)
    end

    ## needle areas
    total_needle_area = 50.8 * 1e2^2 # cm^2; from m^2 (input)
    needle_areas = [i * total_needle_area / sum(1:n_branches)
					for i in n_branches:-1:1]

    ## branchgrowth frequency
    branchgrowth_frequency = 3 
		# how many steps it takes branches to grow a new segment (estimate)

    ## root dimensions
    total_root_area = 96.1 * 1e2^2 # cm^2; from m^2 (input)

    # parameters

    ## hydraulic conductivities
    permeability = 2.5e-12 * 1e2^2 # cm^2; from m^2 (input)
    η = 1.0 * 1e-9 * (1/3600) # MPa hr; from mPa s (dynamic viscosity of water)
    l = segment_length # cm
    ρ_w = 1.0 # g / cm^3 (density of water)
    K_s_stem = ρ_w * permeability / η / l # g / hr / MPa / cm^2

    L_p = 3.2e-8 * 1e2 * 3600 # cm / hr / MPa; from m / s / MPa (input)
    K_roots = ρ_w * L_p * total_root_area # g / hr / MPa

    ## elastic modulus
    E_D_stem = [0.15 * 1e3, 17.5 * 0.15 * 1e3] # MPa; from GPa (input)
end

struct HollowCylinder{T} <: ModuleShape
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

mutable struct Stem <: Node
    D::Vector
    h
end

mutable struct Branch <: Node
    D::Vector
    h
end

mutable struct StemTip <: Node
    D::Vector
    h
    branch_nr
    age
end

mutable struct BranchTip <: Node
    D::Vector
    h
    age
    needle_area
end

# ## Instantiation


FINAL_STEM_RX04 = 2.0
plant = +([Stem([sapwood_r0, segment_length], h) for h in 0.0:0.2:FINAL_STEM_RX04]...) + 
    Branch([sapwood_r0/2, segment_length], FINAL_STEM_RX04 + 0.2) +
    BranchTip([sapwood_r0/4, segment_length], FINAL_STEM_RX04 + 0.2, 0, needle_areas[1])

struct Soil <: Node end
struct Air <: Node end

graphs = [Soil(), plant, Air()]
intergraph_connections = [(1, 2) => (:Soil, getnodes(plant)[1]), (2, 3) => (:BranchTip, :Air)]
plantstructure = PlantStructure(graphs, intergraph_connections)

# # Function

# ## New functional modules

import PlantModules: t, d
readcsv(filepath) = readlines(filepath) .|> x -> split(x, ", ") .|> x -> parse(Float64, x)

# ### needle area

function needle_area_module(; name, needle_area)
    @parameters begin
        needle_area_val = needle_area, [description = "Value of total needle area on branch"]
    end 
    @variables begin
        needle_area(t), [description = "Total needle area on branch"]
    end
    eqs = [needle_area ~ needle_area_val]
    return System(eqs, t; name)
end

# ### Fixed transpiration

# #### transpiration rate from data
transpiration_data = readcsv("./tutorials/clouddata/transpiration_data.csv")
transpiration_rate = LinearInterpolation(last.(transpiration_data), first.(transpiration_data))

# #### module
function fixed_daynight_hydraulic_connection(; name, original_order)
    @constants begin
        t_unit = 1, [description = "Dummy constant for correcting units"]
        uc = 1e-3 * 3600 * (1e-2)^2, [description = "Unit conversion from (mg / s / m^2) to (g / hr / cm^2)"]
    end
    @variables begin
        F_s(t), [description = "Specific water flux from compartment 2 to compartment 1"]
        F(t), [description = "Water flux from compartment 2 to compartment 1"]
        needle_area(t), [description = "Total needle area on branch"]
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
	:Stem => [hydraulic_module, constant_carbon_module, K_module],
    :Branch => [hydraulic_module, constant_carbon_module, K_module],
	:BranchTip => [hydraulic_module, needle_area_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module],
)

connecting_modules = Dict(
	(:Soil, :Stem) => constant_hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Branch) => hydraulic_connection,
    (:Branch, :Branch) => hydraulic_connection,
	(:Branch, :BranchTip) => hydraulic_connection,
	(:BranchTip, :Air) => fixed_daynight_hydraulic_connection,
)

plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# # Parameters
begin
    default_changes = Dict(
		:needle_area => 0.0,
		:E_D => 1.0,#E_D_stem,
        :K_s => K_s_stem,
        :ϕ_D => 1.0, :M => 0.0, 
		:Ψ => PlantModules.soilfunc(0.9), 
		:shape => HollowCylinder(radial_fraction_sapwood)
	)

	module_defaults = Dict(
        :Soil => Dict(:W_max => 1e6, :W_r => 0.9),
    )

    connection_values = Dict(
        (:Soil, :Stem) => Dict(:K => K_roots),
    )

    plantparams = PlantParameters(; default_changes, module_defaults, connection_values)
end;

# # Run it

tspan = (0.0, 24.0);
system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false);
@time prob = ODEProblem(system, [], tspan, sparse = true, use_scc = false);
@time sol = solve(prob, FBDF());

using Plots
plotgraph(sol, plantstructure, structmod = :Stem, varname = :W)

# # Show it

## local sensitivity

stem_radius_var = get_subsystem_variables(system, plantstructure, :D, :Stem)[1][2]
K_s_vars = get_subsystem_variables(system, plantstructure, :K_s, [:Stem, :Branch])
function get_diameter(K_s)
    newprob = remake(prob, p = Pair.(K_s_vars, [K_s]))
    # newprob = remake(newprob, 
    #     p = MTKParameters(typeof(K_s).(newprob.p[1]), typeof(K_s).(newprob.p[2]), (), (typeof(K_s).(newprob.p[3]),), (), ()),
    #     u0 = typeof(K_s).(newprob.u0)
    # )
	newsol = solve(newprob, FBDF(), saveat = 0.1)
    return newsol[stem_radius_var]
end

using ForwardDiff

Base.Float64(x::ForwardDiff.Dual) = x.value #!

K_sensitivity = ForwardDiff.derivative(get_diameter, K_s_stem)

plot(K_sensitivity)