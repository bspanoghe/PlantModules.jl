### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 34c2f228-076d-4d1b-ac90-6d9bf683628b
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("../..")

# ╔═╡ e9c09211-04b6-464b-9d0b-3a33bef74073
using PlutoUI; TableOfContents()

# ╔═╡ 792edd4e-aec4-48d6-801f-8a48de61e11b
using PlantModules

# ╔═╡ f6ebbd49-68d2-4ae6-92c3-aebe38e5f4ac
using ModelingToolkit, OrdinaryDiffEq, Plots

# ╔═╡ 1966353a-462a-4784-91fa-8a6929dc81a9
using PlantGraphs

# ╔═╡ 3c306ef3-2546-435e-b1b6-c5325499590e
md"# Package speed benchmarking"

# ╔═╡ 1b98b39f-c3b9-49f2-9f34-ee854b2da83c
md"""
In this notebook, we will benchmark the speed of the package. Considering our framework can technically create any DAE-based model, we will limit this benchmark to the core functional modules of `PlantModules.jl`. Additionally, we consider two very straightforward plant structures for simplicity's sake.
"""

# ╔═╡ 02479490-ddee-4f5c-b34a-c6225c449a9e
md"## Setup"

# ╔═╡ 92d46933-9ba5-4e10-89ef-2cfa1e88ed6a
md"## System definition"

# ╔═╡ bb2dc233-c5cc-49c6-aef1-d46e2f0cff0c
md"### Structure"

# ╔═╡ 7e67ea9b-2e78-462f-93a4-5bbc3a50e1e8
md"""
For plant structure, we will consider plants of a single node type `Segment` connected to an environment consisting of a single soil and air node. In order to inspect the influence of the level of branching on solution time, we consider two types of structure:
- A completely linear plant, where the base is connected to the soil and all other nodes are connected to the air. We can imagine this structure to represent something like a conifer, with each plant segment containing a part of the stem and all connected branches.
- A plant of which each segment branches into two other segments, where the base is connected to the soil and all leaf nodes are connected to the air. We can imagine this to represent a herbaceous plant of some kind.
"""

# ╔═╡ ca93a155-e8cd-4fa0-91cb-7c9d5c960c74
md"""
We will define the structure of the plants using rewrite rules so we can easily change the size of the system, using `PlantGraphs.jl` from the [`VirtualPlantLab.jl`](https://virtualplantlab.com/stable/) ecosystem.
"""

# ╔═╡ 67ba4b28-46bd-4863-aeaf-e57ac456068f
Base.@kwdef struct Segment <: PlantGraphs.Node
	D::Vector{Float64} = PlantModules.default_values[:D]
end

# ╔═╡ ab8eb6b5-c297-4f23-b7b9-681c7a5a8a75
struct Soil <: PlantGraphs.Node end

# ╔═╡ 68f161ea-ff5a-46df-8b65-929f27265c31
struct Air <: PlantGraphs.Node end

# ╔═╡ 10098f8c-fdb2-4ebc-85e7-f9212874ac74
md"#### Linear"

# ╔═╡ 038c3d32-0c12-495d-81c8-175c3b13a038
# all segments split into two linearly connected segments
linear_rule = Rule(
	Segment,
	rhs = seg -> Segment() + Segment()
)

# ╔═╡ 3d18d317-ff34-4bcd-83d9-479db8305005
# get the structure of the linear plant (without and with connections to environment) for a given number of rewrite steps
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

# ╔═╡ 9e269dd6-4eb0-4812-9349-b1e203f92cd3
md"#### Branching"

# ╔═╡ d434d804-6f6c-442b-8130-c61f6994701b
# all leaf nodes branch out into two more segments, where the cross area of the parent segment equals the summed cross areas of the child segments (as per Da Vinci's rule)
branching_rule = Rule(
	Segment,
	lhs = seg -> !has_children(seg),
	rhs = seg -> Segment(data(seg).D) + 
		(
			Segment(data(seg).D .* [1/sqrt(2), 1]), 
			Segment(data(seg).D .* [1/sqrt(2), 1])
		)
)

# ╔═╡ 73d72fb9-82b9-44c1-b56c-4bd7475bafb6
# get the structure of the branching plant (without and with connections to environment) for a given number of rewrite steps
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

# ╔═╡ cb1628d4-d85d-4f5c-b9ad-2ae29f3df899
md"### Coupling"

# ╔═╡ 4260b406-8d16-460a-9d24-76c66f68fa4f
md"""
All node types are assigned default hydraulic functionality.
"""

# ╔═╡ e910c522-4aaf-403d-9f2d-37d7d011046c
module_coupling = Dict(
	:Segment => [hydraulic_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
);

# ╔═╡ 51958737-a4f1-43af-98f5-7077418345df
connecting_modules = Dict(
	(:Soil, :Segment) => constant_hydraulic_connection,
	(:Segment, :Segment) => hydraulic_connection,
	(:Segment, :Air) => daynight_hydraulic_connection
);

# ╔═╡ 74ad2d32-27b4-4c99-8132-37d208869a98
const plantcoupling = PlantCoupling(; module_coupling, connecting_modules); # mark global variables as having a constant type for increased efficiency when using them in functions

# ╔═╡ b1b8a817-b2e5-4aa3-93f1-62dcc48a8787
md"### Parameters"

# ╔═╡ 0f946322-28ac-4924-bb8f-f0c7eebb4c15
md"""
The parameter values are defined in function of the size of the plant in order to scale the amount of water in the soil to the amount of transpiration taking place. This is to ensure all simulations have a similar profile for soil water potential through time, as this can have a significant effect on solving time. All other parameters are set to their default values, with two exceptions:
- The specific hydraulic conductivity `K_s` of the plant segments is increased to prevent larger systems from having unrealistically small hydraulic conductivities, considering the cross area of their segments is not increased in our model.
- The hydraulic conductivity `K` of the air is set to a smaller, more realistic value, as the default value is intended to represent a well-conducting node.
- The relative water content `W_r` of the air is set to a smaller value for similar reasons.
"""

# ╔═╡ 51003cf2-0d80-4a40-bb52-cfec7920c798
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

# ╔═╡ 3fea0aca-d6c2-43b8-b852-25bf8276f33a
md"## Sanity check"

# ╔═╡ bb0a5076-6f08-4eab-a762-e9231ad03949
md"""
Before getting to the actual benchmarks, we'll do a quick sanity check to verify our system definition returns sufficiently realistic simulations.
"""

# ╔═╡ 4cf1bd36-5742-41f6-916f-c7a173793920
function test_simulation(plantstructure; tspan = (0.0, 7*24.0))
	plantparams = get_params(plantstructure)
	system = generate_system(plantstructure, plantcoupling, plantparams)
	prob = ODEProblem(system, [], tspan, sparse = true)
	sol = solve(prob, FBDF())

	return sol
end

# ╔═╡ d0b30517-7eee-4573-8081-58b754842c1e
rewrite_steps_test = 5

# ╔═╡ 959532f7-12c4-4af8-9e33-5a365a468966
md"### Linear"

# ╔═╡ 631a364f-7e27-4310-9149-9f0938191dc2
test_plant_linear, test_structure_linear = get_structure_linear(rewrite_steps_test)

# ╔═╡ 1fae811b-32a7-451a-9eae-d872da8ace68
Plots.plot(
	plotstructure(test_plant_linear, title = "Without environment"),
	plotstructure(test_structure_linear, title = "With environment"),
	margins = 5*Plots.mm, size = (800, 400)
)

# ╔═╡ c23b0d2c-a051-4303-80e3-1082ea9a94bc
test_sol_linear = test_simulation(test_structure_linear);

# ╔═╡ 186b6c66-85ec-4260-9e11-7370e8773a0b
begin
	plotgraph(test_sol_linear, test_structure_linear, 
			  varname = :P, structmod = :Segment, 
			  label = "Segment turgor pressure", 
			  title = "Turgor pressure over time",
			  xlabel = "Time (h)", ylabel = "P (MPa)",
			  size = (800, 400), margins = 5*Plots.mm,
			  linewidth = 2
	)
	hline!(
		[PlantModules.default_values[:Γ]], label = "Yield turgor pressure", lw = 2
	)
end

# ╔═╡ 0d1fbc8b-dbae-4a25-8806-0f917ea21f33
plotgraph(
	test_sol_linear, test_structure_linear, varname = :W, structmod = :Soil,
	title = "Soil water content over time", legend = false,
	xlabel = "Time (h)", ylabel = "Water content (g)",
	size = (800, 400), margins = 5*Plots.mm, lw = 2
)

# ╔═╡ 16de2ad2-96d4-42f3-a7f1-981059dc2611
md"### Branching"

# ╔═╡ 53d77b90-49cb-44c0-8ef6-24157fa5abe8
test_plant_branching, test_structure_branching = get_structure_branching(rewrite_steps_test)

# ╔═╡ 8fb7f517-e071-4c0d-b739-a4825708426a
Plots.plot(
	plotstructure(test_plant_branching, title = "Without environment"),
	plotstructure(test_structure_branching, title = "With environment"),
	margins = 5*Plots.mm, size = (800, 400)
)

# ╔═╡ cea4d360-f9f6-4b48-96fd-7c2cc82a0d77
test_sol_branching = test_simulation(test_structure_branching);

# ╔═╡ 01be5dc8-8c1f-497b-8f8f-a92ac3c5b29d
begin
	plotgraph(test_sol_branching, test_structure_branching, 
			  varname = :P, structmod = :Segment, 
			  label = "Segment turgor pressure", 
			  title = "Turgor pressure over time",
			  xlabel = "Time (h)", ylabel = "P (MPa)",
			  size = (800, 400), margins = 5*Plots.mm,
			  linewidth = 2
	)
	hline!(
		[PlantModules.default_values[:Γ]], label = "Yield turgor pressure", lw = 2
	)
end

# ╔═╡ 61e5f480-bfc8-470d-bc0f-50853b2c08ca
plotgraph(
	test_sol_branching, test_structure_branching, varname = :W, structmod = :Soil,
	title = "Soil water content over time", legend = false,
	xlabel = "Time (h)", ylabel = "Water content (g)",
	size = (800, 400), margins = 5*Plots.mm, lw = 2
)

# ╔═╡ 4d43634f-b03a-4988-997e-d133e74dca2c
md"## Benchmarking"

# ╔═╡ db53873c-5213-4b66-bcf1-0213b522cab8
md"""
Now let's get to the actual benchmarking. We will time the system generation, problem generation and problem solving for differing numbers of rewrite steps to evaluate how computation time scales with system size.
"""

# ╔═╡ d491c24e-46e8-40ca-9f5f-2f76e7e95c06
function get_stats(plantstructure; tspan = (0.0, 7*24.0))
	plantparams = get_params(plantstructure)
	
	system_stats = @timed generate_system(
		plantstructure, plantcoupling, plantparams
	)
	prob_stats = @timed ODEProblem(
		system_stats.value, [], tspan, sparse = true
	)
	sol_stats = @timed solve(prob_stats.value, FBDF())

	num_nodes = length(getnodes(plantstructure))
	num_variables = length(unknowns(system_stats.value))

	return num_nodes, num_variables, system_stats, prob_stats, sol_stats
end

# ╔═╡ 49b41f29-0b19-483d-b903-70a254b1f19e
md"""
!!! warning
	Note that all  benchmarks are based on only a single run, which makes them somewhat inconsistent compared to taking the median of a sample of runs. The reason for this is that we also benchmark function compilation time, which only triggers the first time a function method is used, making it non-trivial to get a sample of.
"""

# ╔═╡ 45de608d-9bd1-49c8-8dd3-1d97b47ad2bc
function benchmark_linear(rewrite_steps)
	plant, plantstructure = get_structure_linear(rewrite_steps)
	return get_stats(plantstructure)
end

# ╔═╡ c484bbd7-6171-4976-8d0a-a92e2bf7980a
function benchmark_branching(rewrite_steps)
	plant, plantstructure = get_structure_branching(rewrite_steps)
	return get_stats(plantstructure)
end

# ╔═╡ 9900c521-75e0-492c-9010-24034c6ddf62
function stackedbar(xs, ys; kwargs...)
	same_x_groups = [findall(x -> x == u, xs) for u in unique(xs)]
	xs_grouped = [xs[group] for group in same_x_groups]
	ys_grouped = [ys[group] for group in same_x_groups]
	ys_grouped_cs = [cumsum(ys_group) for ys_group in ys_grouped]
	ys_grouped_lower = [[0; ys_group_cs[1:end-1]] for ys_group_cs in ys_grouped_cs]

	xs_flat = reduce(vcat, xs_grouped) |> 
		x -> reshape(x, :, length(unique(xs))) |>
			permutedims
	ys_cs_flat = reduce(vcat, ys_grouped_cs) |> 
		x -> reshape(x, :, length(unique(xs))) |>
			permutedims
	ys_lower_flat = reduce(vcat, ys_grouped_lower) |>
		x -> reshape(x, :, length(unique(xs))) |>
			permutedims
	
	bar(xs_flat, ys_cs_flat, fillto = ys_lower_flat; kwargs...)
end

# ╔═╡ d7a782ee-1d2b-4101-9ad6-dc8747271c49
function plot_size(rewrite_steps_set, stats)
	Plots.plot(
		stackedbar(rewrite_steps_set, getindex.(stats, 1), label = false, 
				   xlabel = "Rewrite steps", ylabel = "Number of nodes",
				   xticks = rewrite_steps_set),
		stackedbar(rewrite_steps_set, getindex.(stats, 2), label = false, 
				   xlabel = "Rewrite steps", ylabel = "Number of variables", 
				   xticks = rewrite_steps_set),
		size = (1000, 400), margins = 5*Plots.mm, plot_title = "System size in function of rewrite steps"
	)
end

# ╔═╡ e191e8cd-94fe-4b90-9a99-b31c6a6251f1
function plot_time(rewrite_steps_set, stats)
	Plots.plot(
		stackedbar(
			reduce(vcat, [fill(i, 3) for i in rewrite_steps_set]),
			reduce(vcat, [[stat.time for stat in stats[i][3:5]] for i in eachindex(stats)]),
			label = ["System generation" "Problem generation" "Solving"],
			xlabel = "Rewrite steps", ylabel = "Computation time (s)",
			title = "Including compilation", xticks = rewrite_steps_set
		),
		stackedbar(
			reduce(vcat, [fill(i, 3) for i in rewrite_steps_set]),
			reduce(vcat, [[stat.time - stat.compile_time for stat in stats[i][3:5]] for i in eachindex(stats)]),
			label = ["System generation" "Problem generation" "Solving"],
			xlabel = "Rewrite steps", ylabel = "Computation time (s)",
			title = "Excluding compilation", xticks = rewrite_steps_set
		),
		size = (1000, 400), margins = 5*Plots.mm, plot_title = "Computation time in function of rewrite steps", plot_titlevspan = 0.1
	)
end

# ╔═╡ 52db86d9-fdff-43b7-978a-19553d218afc
rewrite_steps_set = [2, 4, 6] # note not to include `rewrite_steps_test`, as the functions have already compiled for this system size

# ╔═╡ f4d4311a-066e-4111-b134-0e465a41970d
md"### Linear structure"

# ╔═╡ 1e339ee4-7b19-4bb2-b2f1-7f82d9cd1762
stats_linear = [benchmark_linear(rewrite_steps) for rewrite_steps in rewrite_steps_set];

# ╔═╡ 2a3be0b8-3b02-49d6-ab14-b6c5482891a2
p_size_linear = plot_size(rewrite_steps_set, stats_linear)

# ╔═╡ 2499e650-b54c-4158-9339-21f846f2695e
p_time_linear = plot_time(rewrite_steps_set, stats_linear)

# ╔═╡ 27431fdd-8159-4048-aaf7-084489cd9d29
md"### Branching structure"

# ╔═╡ 7d3047b5-c8ff-4c07-8b27-e76b18534c29
stats_branching = [benchmark_branching(rewrite_steps) for rewrite_steps in rewrite_steps_set];

# ╔═╡ 347160b8-368a-420a-8f65-8b50e4449c2d
p_size_branching = plot_size(rewrite_steps_set, stats_branching)

# ╔═╡ c0cedaa5-5ca2-4479-bead-7ee289cc5fcf
p_time_branching = plot_time(rewrite_steps_set, stats_branching)

# ╔═╡ Cell order:
# ╟─3c306ef3-2546-435e-b1b6-c5325499590e
# ╟─1b98b39f-c3b9-49f2-9f34-ee854b2da83c
# ╟─02479490-ddee-4f5c-b34a-c6225c449a9e
# ╠═34c2f228-076d-4d1b-ac90-6d9bf683628b
# ╠═e9c09211-04b6-464b-9d0b-3a33bef74073
# ╠═792edd4e-aec4-48d6-801f-8a48de61e11b
# ╠═f6ebbd49-68d2-4ae6-92c3-aebe38e5f4ac
# ╠═1966353a-462a-4784-91fa-8a6929dc81a9
# ╟─92d46933-9ba5-4e10-89ef-2cfa1e88ed6a
# ╟─bb2dc233-c5cc-49c6-aef1-d46e2f0cff0c
# ╟─7e67ea9b-2e78-462f-93a4-5bbc3a50e1e8
# ╟─ca93a155-e8cd-4fa0-91cb-7c9d5c960c74
# ╠═67ba4b28-46bd-4863-aeaf-e57ac456068f
# ╠═ab8eb6b5-c297-4f23-b7b9-681c7a5a8a75
# ╠═68f161ea-ff5a-46df-8b65-929f27265c31
# ╟─10098f8c-fdb2-4ebc-85e7-f9212874ac74
# ╠═038c3d32-0c12-495d-81c8-175c3b13a038
# ╠═3d18d317-ff34-4bcd-83d9-479db8305005
# ╟─9e269dd6-4eb0-4812-9349-b1e203f92cd3
# ╠═d434d804-6f6c-442b-8130-c61f6994701b
# ╠═73d72fb9-82b9-44c1-b56c-4bd7475bafb6
# ╟─cb1628d4-d85d-4f5c-b9ad-2ae29f3df899
# ╟─4260b406-8d16-460a-9d24-76c66f68fa4f
# ╠═e910c522-4aaf-403d-9f2d-37d7d011046c
# ╠═51958737-a4f1-43af-98f5-7077418345df
# ╠═74ad2d32-27b4-4c99-8132-37d208869a98
# ╟─b1b8a817-b2e5-4aa3-93f1-62dcc48a8787
# ╟─0f946322-28ac-4924-bb8f-f0c7eebb4c15
# ╠═51003cf2-0d80-4a40-bb52-cfec7920c798
# ╟─3fea0aca-d6c2-43b8-b852-25bf8276f33a
# ╟─bb0a5076-6f08-4eab-a762-e9231ad03949
# ╠═4cf1bd36-5742-41f6-916f-c7a173793920
# ╠═d0b30517-7eee-4573-8081-58b754842c1e
# ╟─959532f7-12c4-4af8-9e33-5a365a468966
# ╠═631a364f-7e27-4310-9149-9f0938191dc2
# ╟─1fae811b-32a7-451a-9eae-d872da8ace68
# ╠═c23b0d2c-a051-4303-80e3-1082ea9a94bc
# ╟─186b6c66-85ec-4260-9e11-7370e8773a0b
# ╟─0d1fbc8b-dbae-4a25-8806-0f917ea21f33
# ╟─16de2ad2-96d4-42f3-a7f1-981059dc2611
# ╠═53d77b90-49cb-44c0-8ef6-24157fa5abe8
# ╟─8fb7f517-e071-4c0d-b739-a4825708426a
# ╠═cea4d360-f9f6-4b48-96fd-7c2cc82a0d77
# ╟─01be5dc8-8c1f-497b-8f8f-a92ac3c5b29d
# ╟─61e5f480-bfc8-470d-bc0f-50853b2c08ca
# ╟─4d43634f-b03a-4988-997e-d133e74dca2c
# ╟─db53873c-5213-4b66-bcf1-0213b522cab8
# ╠═d491c24e-46e8-40ca-9f5f-2f76e7e95c06
# ╟─49b41f29-0b19-483d-b903-70a254b1f19e
# ╟─45de608d-9bd1-49c8-8dd3-1d97b47ad2bc
# ╟─c484bbd7-6171-4976-8d0a-a92e2bf7980a
# ╟─9900c521-75e0-492c-9010-24034c6ddf62
# ╟─d7a782ee-1d2b-4101-9ad6-dc8747271c49
# ╟─e191e8cd-94fe-4b90-9a99-b31c6a6251f1
# ╠═52db86d9-fdff-43b7-978a-19553d218afc
# ╟─f4d4311a-066e-4111-b134-0e465a41970d
# ╠═1e339ee4-7b19-4bb2-b2f1-7f82d9cd1762
# ╟─2a3be0b8-3b02-49d6-ab14-b6c5482891a2
# ╟─2499e650-b54c-4158-9339-21f846f2695e
# ╟─27431fdd-8159-4048-aaf7-084489cd9d29
# ╠═7d3047b5-c8ff-4c07-8b27-e76b18534c29
# ╟─347160b8-368a-420a-8f65-8b50e4449c2d
# ╟─c0cedaa5-5ca2-4479-bead-7ee289cc5fcf
