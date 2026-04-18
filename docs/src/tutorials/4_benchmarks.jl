### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 34c2f228-076d-4d1b-ac90-6d9bf683628b
using Pkg; Pkg.activate("../..")

# ╔═╡ e9c09211-04b6-464b-9d0b-3a33bef74073
using PlutoUI; TableOfContents()

# ╔═╡ 792edd4e-aec4-48d6-801f-8a48de61e11b
using PlantModules

# ╔═╡ f6ebbd49-68d2-4ae6-92c3-aebe38e5f4ac
using ModelingToolkit, OrdinaryDiffEq, Plots

# ╔═╡ 1966353a-462a-4784-91fa-8a6929dc81a9
using PlantGraphs

# ╔═╡ 92d46933-9ba5-4e10-89ef-2cfa1e88ed6a
md"## System setup"

# ╔═╡ bb2dc233-c5cc-49c6-aef1-d46e2f0cff0c
md"### Structural modules"

# ╔═╡ 67ba4b28-46bd-4863-aeaf-e57ac456068f
struct Segment <: PlantGraphs.Node end

# ╔═╡ ab8eb6b5-c297-4f23-b7b9-681c7a5a8a75
struct Soil <: PlantGraphs.Node end

# ╔═╡ 68f161ea-ff5a-46df-8b65-929f27265c31
struct Air <: PlantGraphs.Node end

# ╔═╡ 17f47e5c-1040-4252-85bf-5a8e5dba4119
md"### Structural connections"

# ╔═╡ 10098f8c-fdb2-4ebc-85e7-f9212874ac74
md"#### Linear"

# ╔═╡ 038c3d32-0c12-495d-81c8-175c3b13a038
linear_rule = Rule(
	Segment,
	rhs = seg -> Segment() + Segment()
)

# ╔═╡ 3d18d317-ff34-4bcd-83d9-479db8305005
function get_structure_linear(rewrite_steps)
	plant = Graph(axiom = Segment(), rules = linear_rule)
    for i in 1:rewrite_steps
        rewrite!(plant)
    end

	plantstructure = PlantStructure(
		[plant, Soil(), Air()],
		[
			(1, 2) => (getnodes(plant)[1], :Soil),
			(1, 3) => (seg, air) -> seg != getnodes(plant)[1]
		]
	)
	
    return plant, plantstructure
end

# ╔═╡ 9e269dd6-4eb0-4812-9349-b1e203f92cd3
md"#### Branching"

# ╔═╡ d434d804-6f6c-442b-8130-c61f6994701b
branching_rule = Rule(
	Segment,
	lhs = seg -> !has_children(seg),
	rhs = seg -> Segment() + (Segment(), Segment()) #! da vinci
)

# ╔═╡ 73d72fb9-82b9-44c1-b56c-4bd7475bafb6
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
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ 3fea0aca-d6c2-43b8-b852-25bf8276f33a
md"## Sanity check"

# ╔═╡ 0e302e9f-84b5-40c4-8d3d-193c92a0947e
@bind get_structure Select([get_structure_linear, get_structure_branching])

# ╔═╡ a4b12858-b40d-41a8-82a8-2e69b91d94ab
plant, plantstructure = get_structure(3)

# ╔═╡ 1fae811b-32a7-451a-9eae-d872da8ace68
Plots.plot(
	plotstructure(plant),
	plotstructure(plantstructure),
	margins = 5*Plots.mm
)

# ╔═╡ 8b34a2a9-971c-4dea-9d56-95703defe965
num_evaporating_segments = length(getneighbors(getnodes(plantstructure)[end], plantstructure))

# ╔═╡ 35383f8d-c047-458b-ab3f-cdaf3be3eea9
module_defaults = Dict(
	:Soil => Dict(:W_max => num_evaporating_segments * 1e2, :K => 1.0),
	:Air => Dict(:W_r => 0.7, :K => 1e-2)
);

# ╔═╡ 580d3ebb-2a5c-4762-ab55-d7e271b899c0
plantparams = PlantParameters(; module_defaults);

# ╔═╡ 8589bed6-d602-4ff3-bedb-6d57d8395ed5
system = generate_system(plantstructure, plantcoupling, plantparams);

# ╔═╡ a8d4a450-ca20-48d4-b3e6-b51424a1ffed
tspan = (0.0, 7*24.0)

# ╔═╡ 9ec62911-e0c2-4844-8df3-c15350f20752
prob = ODEProblem(system, [], tspan, sparse = true);

# ╔═╡ 385ce3af-d10e-4cf1-8280-7d64e8e8e0f5
sol = solve(prob, FBDF());

# ╔═╡ 61268c0d-4fcc-4b82-ae9a-d1979d9a3159
plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Soil, :Segment])

# ╔═╡ 186b6c66-85ec-4260-9e11-7370e8773a0b
plotgraph(sol, plantstructure, varname = :W, structmod = [:Soil, :Segment])

# ╔═╡ 4d43634f-b03a-4988-997e-d133e74dca2c
md"## Benchmarking"

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

# ╔═╡ d491c24e-46e8-40ca-9f5f-2f76e7e95c06
function get_stats(plantstructure)	
	num_evaporating_segments = getneighbors(
		getnodes(plantstructure)[end],
		plantstructure
	) |> length
	module_defaults = Dict(
		:Soil => Dict(:W_max => num_evaporating_segments * 1e2, :K => 1.0),
		:Air => Dict(:W_r => 0.7, :K => 1e-2)
	)
	plantparams = PlantParameters(; module_defaults)
	
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

# ╔═╡ 52db86d9-fdff-43b7-978a-19553d218afc
rewrite_steps_set = 1:5

# ╔═╡ d7a782ee-1d2b-4101-9ad6-dc8747271c49
function plot_size(stats)
	Plots.plot(
		stackedbar(rewrite_steps_set, getindex.(stats, 1), label = false, 
				   xlabel = "Rewrite steps", ylabel = "Number of nodes"),
		stackedbar(rewrite_steps_set, getindex.(stats, 2), label = false, 
				   xlabel = "Rewrite steps", ylabel = "Number of variables"),
		size = (1000, 400), margins = 5*Plots.mm, plot_title = "System size in function of rewrite steps"
	)
end

# ╔═╡ e191e8cd-94fe-4b90-9a99-b31c6a6251f1
function plot_time(stats)
	Plots.plot(
		stackedbar(
			reduce(vcat, [fill(i, 3) for i in rewrite_steps_set]),
			reduce(vcat, [[stat.time for stat in stats[i][3:5]] for i in eachindex(stats)]),
			label = ["System generation" "Problem generation" "Solving"],
			xlabel = "Rewrite steps", ylabel = "Computation time (s)",
			title = "Including compilation"
		),
		stackedbar(
			reduce(vcat, [fill(i, 3) for i in rewrite_steps_set]),
			reduce(vcat, [[stat.time - stat.compile_time for stat in stats[i][3:5]] for i in eachindex(stats)]),
			label = ["System generation" "Problem generation" "Solving"],
			xlabel = "Rewrite steps", ylabel = "Computation time (s)",
			title = "Excluding compilation"
		),
		size = (1000, 400), margins = 5*Plots.mm, plot_title = "Computation times in function of rewrite steps", plot_titlevspan = 0.1
	)
end

# ╔═╡ 1e339ee4-7b19-4bb2-b2f1-7f82d9cd1762
stats_linear = [benchmark_linear(rewrite_steps) for rewrite_steps in rewrite_steps_set];

# ╔═╡ 2a3be0b8-3b02-49d6-ab14-b6c5482891a2
plot_size(stats_linear)

# ╔═╡ 2499e650-b54c-4158-9339-21f846f2695e
plot_time(stats_linear)

# ╔═╡ 7d3047b5-c8ff-4c07-8b27-e76b18534c29
stats_branching = [benchmark_branching(rewrite_steps) for rewrite_steps in rewrite_steps_set];

# ╔═╡ 347160b8-368a-420a-8f65-8b50e4449c2d
plot_size(stats_branching)

# ╔═╡ c0cedaa5-5ca2-4479-bead-7ee289cc5fcf
plot_time(stats_branching)

# ╔═╡ Cell order:
# ╠═34c2f228-076d-4d1b-ac90-6d9bf683628b
# ╠═e9c09211-04b6-464b-9d0b-3a33bef74073
# ╠═792edd4e-aec4-48d6-801f-8a48de61e11b
# ╠═f6ebbd49-68d2-4ae6-92c3-aebe38e5f4ac
# ╠═1966353a-462a-4784-91fa-8a6929dc81a9
# ╟─92d46933-9ba5-4e10-89ef-2cfa1e88ed6a
# ╟─bb2dc233-c5cc-49c6-aef1-d46e2f0cff0c
# ╠═67ba4b28-46bd-4863-aeaf-e57ac456068f
# ╠═ab8eb6b5-c297-4f23-b7b9-681c7a5a8a75
# ╠═68f161ea-ff5a-46df-8b65-929f27265c31
# ╟─17f47e5c-1040-4252-85bf-5a8e5dba4119
# ╟─10098f8c-fdb2-4ebc-85e7-f9212874ac74
# ╠═038c3d32-0c12-495d-81c8-175c3b13a038
# ╠═3d18d317-ff34-4bcd-83d9-479db8305005
# ╟─9e269dd6-4eb0-4812-9349-b1e203f92cd3
# ╠═d434d804-6f6c-442b-8130-c61f6994701b
# ╠═73d72fb9-82b9-44c1-b56c-4bd7475bafb6
# ╟─cb1628d4-d85d-4f5c-b9ad-2ae29f3df899
# ╠═e910c522-4aaf-403d-9f2d-37d7d011046c
# ╠═51958737-a4f1-43af-98f5-7077418345df
# ╠═74ad2d32-27b4-4c99-8132-37d208869a98
# ╟─3fea0aca-d6c2-43b8-b852-25bf8276f33a
# ╠═0e302e9f-84b5-40c4-8d3d-193c92a0947e
# ╠═a4b12858-b40d-41a8-82a8-2e69b91d94ab
# ╠═1fae811b-32a7-451a-9eae-d872da8ace68
# ╠═8b34a2a9-971c-4dea-9d56-95703defe965
# ╠═35383f8d-c047-458b-ab3f-cdaf3be3eea9
# ╠═580d3ebb-2a5c-4762-ab55-d7e271b899c0
# ╠═8589bed6-d602-4ff3-bedb-6d57d8395ed5
# ╠═a8d4a450-ca20-48d4-b3e6-b51424a1ffed
# ╠═9ec62911-e0c2-4844-8df3-c15350f20752
# ╠═385ce3af-d10e-4cf1-8280-7d64e8e8e0f5
# ╠═61268c0d-4fcc-4b82-ae9a-d1979d9a3159
# ╠═186b6c66-85ec-4260-9e11-7370e8773a0b
# ╟─4d43634f-b03a-4988-997e-d133e74dca2c
# ╟─9900c521-75e0-492c-9010-24034c6ddf62
# ╠═d491c24e-46e8-40ca-9f5f-2f76e7e95c06
# ╠═45de608d-9bd1-49c8-8dd3-1d97b47ad2bc
# ╠═c484bbd7-6171-4976-8d0a-a92e2bf7980a
# ╠═d7a782ee-1d2b-4101-9ad6-dc8747271c49
# ╠═e191e8cd-94fe-4b90-9a99-b31c6a6251f1
# ╠═52db86d9-fdff-43b7-978a-19553d218afc
# ╠═1e339ee4-7b19-4bb2-b2f1-7f82d9cd1762
# ╠═2a3be0b8-3b02-49d6-ab14-b6c5482891a2
# ╠═2499e650-b54c-4158-9339-21f846f2695e
# ╠═7d3047b5-c8ff-4c07-8b27-e76b18534c29
# ╠═347160b8-368a-420a-8f65-8b50e4449c2d
# ╠═c0cedaa5-5ca2-4479-bead-7ee289cc5fcf
