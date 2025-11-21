using Infiltrator
using Revise
using Plots
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful
using DataInterpolations, ForwardDiff, Sparspak

# for figures
include(homedir() * raw"\Documents\Github\Caverns_of_code\Julia\Lifehacks\catpuccin\get_palette.jl")
cpalette = get_palette("prettycolors")

plotdir = homedir() * raw"\Documents\Github\Den_of_evil\Non-note files\images\\"

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
    ϵ_D_stem = [0.15 * 1e3, 17.5 * 0.15 * 1e3] # MPa; from GPa (input)
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

# ## Growing rules

branching_rule = Rule( # top stem segment grows another thinner stem segment and a branch every `branching_frequency` steps, if stem is old enough to form crown
    StemTip,
    lhs = st -> data(st).age >= crown_start_age && data(st).age % branching_frequency == 0,
    rhs = st -> Stem(data(st).D, data(st).h) + (
        BranchTip([branch_radius(data(st).D[1], Δr), segment_length], data(st).h, 0, needle_areas[data(st).branch_nr]),
        StemTip(data(st).D - [Δr, 0.0], data(st).h + segment_length, data(st).branch_nr + 1, data(st).age + 1)
    )
)

vertical_growth_rule = Rule( # top stem segment grows only another stem segment if not branching
    StemTip,
    lhs = st -> !(data(st).age >= crown_start_age && data(st).age % branching_frequency == 0),
    rhs = st -> Stem(data(st).D, data(st).h) + 
        StemTip(data(st).D, data(st).h + segment_length, data(st).branch_nr, data(st).age + 1)
)

branchgrowth_rule = Rule( # top branch segments grown another branch segment every `branchgrowth_frequency` steps
    BranchTip,
    lhs = br -> data(br).age % branchgrowth_frequency == 0,
    rhs = br -> Branch(data(br).D, data(br).h) +
        BranchTip(data(br).D, data(br).h, data(br).age + 1, data(br).needle_area)
)

branchage_rule = Rule( # branch tips also age if not growing
    BranchTip,
    lhs = br -> !(data(br).age % branchgrowth_frequency == 0),
    rhs = br -> BranchTip(data(br).D, data(br).h, data(br).age + 1, data(br).needle_area)
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

struct Soil <: Node end
struct Air <: Node end

graphs = [Soil(), plant, Air()]
intergraph_connections = [(1, 2) => (:Soil, getnodes(plant)[1]), (2, 3) => (:BranchTip, :Air)]
plantstructure = PlantStructure(graphs, intergraph_connections)
plotstructure(plantstructure, names = "")

begin
    using Random
    Random.seed!(7)
    plot_structure = PlantStructure([Soil(), plant], [(1, 2) => (:Soil, getnodes(plant)[1])])

    plotstructure(plot_structure, 
        names = "", size = (600, 400), axisbuffer = 0.05, label = false, legend = true, 
        palette = cpalette
    )
    plot!([], [], markershape = :hexagon, color = cpalette[1], label = "Soil")
    plot!([], [], markershape = :hexagon, color = cpalette[2], label = "Stem")
    plot!([], [], markershape = :hexagon, color = cpalette[3], label = "Branch")
    plot!([], [], markershape = :hexagon, color = cpalette[4], label = "BranchTip")
    # savefig(plotdir * "fig_plantmodules_ex2_structure.pdf")
end

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

# ### Fixed transpiration

# #### transpiration rate from data
transpiration_data = readcsv("./tutorials/clouddata/transpiration_data.csv")
transpiration_rate = LinearInterpolation(last.(transpiration_data), first.(transpiration_data))
plot(transpiration_rate) # in mg / s / m^2

# #### module
function fixed_daynight_hydraulic_connection(; name, original_order)
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
		:ϵ_D => ϵ_D_stem, :K_s => K_s_stem,
        :ϕ_D => 0.0, :M => 0.0, 
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
end

# # Run it

tspan = (0.0, 24.0)
system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false)
@time prob = ODEProblem(system, [], tspan, sparse = true, use_scc = false)

@time sol = solve(prob, FBDF())

# # Show it

## transpiration
air_water_inflow = get_subsystem_variables(system, plantstructure, :ΣF, :Air)[1]
transpiration_unit_conversion(var) = var * 1e3 / (total_needle_area * (1e-2)^2) / 3600
plot(
    sol, idxs = [transpiration_unit_conversion(air_water_inflow)], label = false,
	xlabel = "Time of day (h)", ylabel = "Transpiration (mg / m^2 / s)",
	xticks = 0:3:24, ylims = (0, 10), yticks = 0:2:10, palette = cpalette
)

## water tension

begin
	stem_pressures = get_subsystem_variables(system, plantstructure, :P, :Stem)
	pressure_segment_heights = [10.0, crown_base_height, 1200.0]
	pressure_segment_nrs = ceil.(Int64, pressure_segment_heights / segment_length)
	
	plot(sol, idxs = stem_pressures[pressure_segment_nrs], ylims = (-0.4, 0.0),
		 xticks = 0.0:3.0:24.0, xlabel = "Time of day (h)", palette = cpalette,
		 ylabel = "Water tension (MPa)",
		 label = ["Base of stem" "Crown base" "Top of tree"])
end

## diameter variation

begin
    dimension_variables = get_subsystem_variables(system, plantstructure, :D, :Stem)

    diameter_change_mm(var, sol) = 2*10*(var - sol[var][1]) 
		# substract first value and change from radius in cm to diameter in mm
	diameter_segment_heights = [6.5 * 1e2, 2.5 * 1e2]
    diameter_segment_nrs = ceil.(Int64, diameter_segment_heights ./ segment_length)
	
    diameter_datas = [readcsv("./tutorials/clouddata/diameter_data_$(idx).csv") 
					  for idx in ("c", "d")]
	
	subplots = [plot(), plot()]
	for i in eachindex(subplots)
		plot!(subplots[i], sol, idxs = 
			  [diameter_change_mm(
				  dimension_variables[diameter_segment_nrs[i]][1], sol
			  )],
			  label = "Simulated", color = cpalette[2])
		plot!(subplots[i], first.(diameter_datas[i]), last.(diameter_datas[i]), 
			  label = "Data", color = cpalette[1])
		plot!(subplots[i], 
			  title = "Stem at height $(Int64(diameter_segment_heights[i])) cm")
	end
	plot(subplots..., layout = (2, 1), ylabel = "Diameter change (mm)", 
		 size = (800, 600), margins = 5*Plots.mm, #ylims = (-0.07, 0.01), 
		 yticks = -0.07:0.01:0.01, xticks = 0:3:24, xlabel = "Time of day (h)")
    # savefig(plotdir * "fig_plantmodules_ex2_results.pdf")
end

plotgraph(sol, plantstructure, structmod = :Branch, varname = :P)

# Uncertainty analysis

D1_at_250cm = dimension_variables[diameter_segment_nrs[2]][1]

## monte carlo

ϵ_D_r_range = [0.03 * 1e3, 0.27 * 1e3]
sample_range(a, b) = (b-a) * rand() + a

begin
    Random.seed!(1337)
	p_montecarlo = plot()
					   
	for _ in 1:20
	    ϵ_D_r_sample = sample_range(ϵ_D_r_range[1], ϵ_D_r_range[2])
	    ϵ_D_sample = [ϵ_D_r_sample, 17.5 * ϵ_D_r_sample] 
			# Perämäki et al. use a constant ratio of 1/17.5 between the radial and longitudinal elastic moduli
	
	    prob2 = remake_graphsystem(
			prob, system, plantstructure, :ϵ_D, 
			[:Stem, :Branch, :BranchTip], ϵ_D_sample
		)
	    @time sol2 = solve(prob2, FBDF(), reltol = 1e-2) 
			# We use a higher relative tolerance for faster solving. Note that this can cause a SingularException error now and again but the faster solving time is generally worth it

	    plot!(p_montecarlo, sol2, 
			  idxs = [diameter_change_mm(D1_at_250cm, sol)], label = false, line_z = ϵ_D_r_sample, color = cgrad(cpalette[[1, 2]]))
	end
	plot!(
		p_montecarlo, xlabel = "Time of day (h)", ylabel = "Diameter change (mm)",
		size = (800, 600), margins = 5*Plots.mm, ylims = (-0.15, 0.01), 
		yticks = -0.15:0.01:0.01, xticks = 0:3:24, title = "Effect of elastic modulus on diameter change"
	)
    # savefig(plotdir * "fig_plantmodules_ex2_montecarlo.pdf")
end

## local sensitivity

function get_diameter(K_s)
    ps = get_subsystem_variables(system, plantstructure, :K_s, [:Stem, :Branch])
    newprob = remake(prob, p = Pair.(ps, [K_s]))
	newsol = solve(newprob, FBDF(), saveat = 0.1)
	segment_diameter = newsol[D1_at_250cm]
    return segment_diameter
end

@time K_sensitivity = ForwardDiff.derivative(get_diameter, K_s_stem)
begin
    plot(
        tspan[1]:0.1:tspan[2], K_sensitivity, 
        xticks = 0:3:24, yticks = 0:0.3e-7:1.5e-7,
        xlabel = "Time of day (h)", ylabel = "Local sensitivity", legend = false, 
        title = "Local sensitivity analysis of the hydraulic conductivity", 
        size = (800, 600), margins = 5*Plots.mm, color = cpalette[1], lw = 2
    )
    # savefig(plotdir * "fig_plantmodules_ex2_sens.pdf")
end