# Prep output file

cd(@__DIR__)

file_name = "benchmark_timings.csv"

if !isfile(file_name)
	open(file_name, "w") do file
		write(file, "type,steps,includes_compilation,sys,prob,sol\n")
	end
end



# Run benchmarks

using Pkg; Pkg.activate("./../docs")
using PlantModules
using ModelingToolkit, OrdinaryDiffEq, Plots
using PlantGraphs

Base.@kwdef struct Segment <: PlantGraphs.Node
	D::Vector{Float64} = PlantModules.default_values[:D]
end

struct Soil <: PlantGraphs.Node end
struct Air <: PlantGraphs.Node end

linear_rule = Rule(
	Segment,
	rhs = seg -> Segment() + Segment()
)

function get_structure_linear(rewrite_steps)
	plant = Graph(axiom = Segment(), rules = linear_rule)
    for i in 1:rewrite_steps
        rewrite!(plant)
    end

	plantstructure = PlantStructure(
		[plant, Soil(), Air()],
		[
			# connect first segment to soil
			(1, 2) => (getnodes(plant)[1], :Soil), 
			# connect all others to the air
			(1, 3) => (seg, air) -> seg != getnodes(plant)[1] 
		]
	)
	
    return plant, plantstructure
end

branching_rule = Rule(
	Segment,
	lhs = seg -> !has_children(seg),
	rhs = seg -> Segment(data(seg).D) + 
		(
			Segment(data(seg).D .* [1/sqrt(2), 1]), 
			Segment(data(seg).D .* [1/sqrt(2), 1])
		)
)

function get_structure_branching(rewrite_steps)
	plant = Graph(axiom = Segment(), rules = branching_rule)
    for i in 1:rewrite_steps
        rewrite!(plant)
    end

	plantstructure = PlantStructure(
		[plant, Soil(), Air()],
		[
			(1, 2) => (getnodes(plant)[1], :Soil),
			(1, 3) => (seg, air) -> is_leaf(seg)
		]
	)
	
    return plant, plantstructure
end

module_coupling = Dict(
	:Segment => [hydraulic_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
);

connecting_modules = Dict(
	(:Soil, :Segment) => constant_hydraulic_connection,
	(:Segment, :Segment) => hydraulic_connection,
	(:Segment, :Air) => daynight_hydraulic_connection
);

const plantcoupling = PlantCoupling(; module_coupling, connecting_modules); # mark global variables as having a constant type for increased efficiency when using them in functions

function get_params(plantstructure; default_changes = Dict{Symbol, Float64}())
	num_evaporating_segments = length(getneighbors(getnodes(plantstructure)[end], plantstructure))
	module_defaults = Dict(
		:Segment => Dict(:K_s => 500.0),
		:Soil => Dict(:W_max => num_evaporating_segments * 1e2),
		:Air => Dict(:K => 1e-2, :W_r => 0.7)
	)
	plantparams = PlantParameters(; default_changes, module_defaults)

	return plantparams
end

function get_stats(plantstructure; tspan = (0.0, 7*24.0))
	plantparams = get_params(plantstructure)
	
	system_stats = @timed generate_system(
		plantstructure, plantcoupling, plantparams
	)
	prob_stats = @timed ODEProblem(
		system_stats.value, [], tspan, sparse = true
	)
	sol_stats = @timed solve(prob_stats.value, FBDF())

	return system_stats, prob_stats, sol_stats
end

function benchmark_linear(rewrite_steps)
	plant, plantstructure = get_structure_linear(rewrite_steps)
	return get_stats(plantstructure)
end

function benchmark_branching(rewrite_steps)
	plant, plantstructure = get_structure_branching(rewrite_steps)
	return get_stats(plantstructure)
end

rewrite_steps_set = [2, 4, 6]

stats_linear = [benchmark_linear(rewrite_steps) for rewrite_steps in rewrite_steps_set];
stats_branching = [benchmark_branching(rewrite_steps) for rewrite_steps in rewrite_steps_set];



# Write results

open(file_name, "a") do file
    for (stats_set, type) in zip([stats_linear, stats_branching], ["linear", "branching"])
        for (stats, steps) in zip(stats_set, rewrite_steps_set)
            for includes_compilation in ["yes", "no"]
                times = [includes_compilation == "yes" ? stat.time : stat.time - stat.compile_time for stat in stats]
                write(file, "$type,$steps,$includes_compilation,$(times[1]),$(times[2]),$(times[3])", "\n")
            end
        end
    end
end