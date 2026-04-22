### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 16e51c70-fe21-40c1-98f0-404254a71b1f
using Pkg

# ╔═╡ ee4fa627-a613-4128-9fb5-5b336a78f18d
begin
	Pkg.activate()
	using Revise
end

# ╔═╡ 6dae54ca-aa7f-4504-b985-195099162109
Pkg.activate("../..")

# ╔═╡ 813b229f-c17e-4a10-9945-dc9ad9066724
using PlutoUI; TableOfContents()

# ╔═╡ b1f03d85-2fa4-4746-a5d8-606671d8375e
using PlantModules

# ╔═╡ 2ea18a37-041a-48a4-9af0-d3e2bc109c63
using ModelingToolkit, OrdinaryDiffEq, Plots

# ╔═╡ 5f494ba3-8c68-4728-9627-3dac9fe7fcd6
using VirtualPlantLab, GLMakie

# ╔═╡ 986648b3-6c8a-465a-87a0-e0ce5fbbefa8
using SkyDomes, PlantBiophysics, 
	PlantBiophysics.PlantMeteo, PlantBiophysics.PlantSimEngine

# ╔═╡ 6b1b3872-8dec-45fa-99d9-9a6903e2170a
using DataInterpolations

# ╔═╡ 2dfdb97f-361c-47bd-b65a-41b02bc8bc57
md"# Tutorial 3: 3D geometry"

# ╔═╡ 87a39d4c-6f22-447e-bd45-1f133b5bb7b3
md"""
In this tutorial, we will more advanced functionality relating to the 3D structure of our plant, including the simulation of different soil compartments and ray tracing. As our package relies on [`VirtualPlantLab.jl`](https://virtualplantlab.com/stable/) for most structural modelling functionality, it will play a key role in this notebook.
"""

# ╔═╡ e52a8761-f703-4dec-b02b-62d3cc831b4f
md"## Setup"

# ╔═╡ dcb943cf-762f-41dd-9231-e07e9d9cda7e
import Plots: plot, plot!

# ╔═╡ 79690228-9815-4bf4-babf-329050d9aca0
import VirtualPlantLab: Mesh

# ╔═╡ d40434c5-4fda-4f2b-a37b-28137924a236
import VirtualPlantLab.PlantGeomPrimitives: Vec

# ╔═╡ a36a2a5b-de0b-46a5-a97c-d0c6b1f2f9f4
md"## Context"

# ╔═╡ 25160d35-aa6b-4ed5-9199-a54d882dbfc9
md"""
As per typical FSPM fashion, we will again model the growth of a tree. However, this time there will be two novelties:
- We simulate carbon dynamics in the leaves based on a ray tracer, a carbon assimilation model and a stomatal conductance model, all of which are available in [`VirtualPlantLab`](https://virtualplantlab.com/stable).
- We simulate water dynamics in a soil discretized into multiple compartments.
"""

# ╔═╡ d3fdc6b7-e69d-425e-ae3b-99e9497b8cb5
md"## Structural definition"

# ╔═╡ 54eb2298-38cf-49dc-8ca9-44872f0d2375
md"""
As we will be using the ray tracing functionality from `VirtualPlantLab.jl`, we recommend users to read the corresponding [tutorial](https://virtualplantlab.com/stable/tutorials/from_tree_forest/raytracedforest/) first, as well as the prior tutorials in the same series. In a nutshell, the extra steps that are required are:
- Define information about the materials of the structural modules that will interact with light;
- Specify this material in the geometry functions;
- Create one or more light sources;
- Add the soil surface to take into account reflected light;
- Create and run the ray tracer.

For didactical purposes, we will base our tutorial directly on the series of tutorials mentioned above. Note that in `VirtualPlantLab.jl` ray tracing tutorial, light interception is used to simulate carbon dynamics based on a discrete sink-source implentation. We will deviate here and instead use `PlantModules.jl`'s functionality to model carbon dynamics coupled to water relations based on differential equations.

The tree structure used in the ray tracing tutorial has a lot of additional information for the discrete sink-source simulation of carbon dynamics. Considering we don't require this information, we will use the simpler tree structure from the [first tutorial](https://virtualplantlab.com/stable/tutorials/from_tree_forest/tree/) in the series.
"""

# ╔═╡ 83fadc9e-74a5-4ddd-915a-0b8692c9280f
md"### Structural modules"

# ╔═╡ 3632551c-2318-4bdf-b813-436c5da7dc68
md"""
The tree consists of the classic set of internodes, nodes and leaves, as well as buds that grow into branches, corresponding nodes left behind by the buds, and meristems that grow into additional phytomers.
"""

# ╔═╡ b27f3786-84b0-469f-ab0c-a60f02030efd
Base.@kwdef struct Meristem <: VirtualPlantLab.Node end

# ╔═╡ b99eb3ee-26dc-4396-9a61-aba45f7462cf
Base.@kwdef struct Bud <: VirtualPlantLab.Node
	D = [0.5] # a small sphere
end

# ╔═╡ 43bb6044-e567-4f20-8be0-5108f5deead8
Base.@kwdef struct Node <: VirtualPlantLab.Node
	D = [0.5, 0.1] # a very short cylinder
end

# ╔═╡ 4ea2dfcf-cc13-46f1-8f70-f42af7a43a09
Base.@kwdef struct BudNode <: VirtualPlantLab.Node
	D = [0.5, 0.1] # a very short cylinder
end

# ╔═╡ fc109086-c42c-474d-8ab0-4add40a07eed
Base.@kwdef mutable struct Internode <: VirtualPlantLab.Node
	D = [0.5, 10] # a long cylinder
	mat::Lambertian{1} = Lambertian(τ = 0.05, ρ = 0.1) # internodes require a material as they interact with light
end

# ╔═╡ 2ff96078-4a91-427b-8c39-689d47600c05
Base.@kwdef mutable struct Leaf <: VirtualPlantLab.Node
	D = [5, 3, 0.05] # a thin cuboid
	mat::Lambertian{1} = Lambertian(τ = 0.05, ρ = 0.1) # leaves also interact with light
	PAR_samples::Vector{Float64} = Float64[]
	PAR_func::Function = zero
end

# ╔═╡ 5a1e5b8a-be42-45e8-9abc-5cecd9dc83ef
md"### Geometry functions"

# ╔═╡ da673a6c-c661-4d00-8e42-365250cb00d3
md"Define a `feed!` method for every structural module that needs a defined geometry, i.e. that needs to be plotted or interacts with light."

# ╔═╡ 1f2b0dc8-84db-4612-841a-6ddaf08345e5
function VirtualPlantLab.feed!(turtle::Turtle, i::Internode, vars)
    # Rotate turtle around the head to implement elliptical phyllotaxis
    rh!(turtle, vars.phyllotaxis)
    HollowCylinder!(
		turtle, length = i.D[2] / 100, height = i.D[1] / 100,
		width = i.D[1] / 100, move = true, colors = RGB(0.5,0.4,0.0), materials = i.mat
	)
    return nothing
end

# ╔═╡ e062618b-3913-4dd8-8ef4-fdcba785fa61
function VirtualPlantLab.feed!(turtle::Turtle, l::Leaf, vars)
    ra!(turtle, -vars.leaf_angle)
    Rectangle!(turtle, length = l.D[1] / 100, width = l.D[2] / 100, move = false,
			 colors = RGB(0.2, 0.6, 0.2), materials = l.mat)
    ra!(turtle, vars.leaf_angle)
	
    return nothing
end

# ╔═╡ d6fa9897-05ae-4388-a299-c08cb62c1414
function VirtualPlantLab.feed!(turtle::Turtle, b::BudNode, vars)
    ra!(turtle, -vars.branch_angle)
end

# ╔═╡ 2c04e08b-1006-4e88-b754-1297eed4dd9c
md"### Structural growth rules"

# ╔═╡ f69abf27-2c51-475b-8344-97361ace66d6
md"""
The tree grows by two rewriting steps:
- Tree meristems grow into phytomers, occuring every rewrite step,
- Buds grow branches, occuring with a probability proportional to the amount of phytomers to the apical meristem,
and one elongation step:
- Every internode elongates for a set fraction of its current length.
"""

# ╔═╡ 120217e8-c708-4a49-8f2c-1191690ce0db
Base.@kwdef struct treeparams # parameters used in rewrite rules
	growth::Float64 = 0.1
	budbreak::Float64 = 0.25
	phyllotaxis::Float64 = 140.0
	leaf_angle::Float64 = 30.0
	branch_angle::Float64 = 45.0
end

# ╔═╡ eb87190d-26e5-4e37-bb09-15b1e89da868
meristem_rule = Rule(
	Meristem,
	rhs = mer -> Node() + (Bud(), Leaf()) + Internode() + Meristem()
)

# ╔═╡ d0e32e1a-733e-4791-82ec-c4a20745a351
function prob_break(bud)
    node = parent(bud)
    check, steps = has_descendant(
		node, 
		condition = n -> data(n) isa Meristem
	)
    steps = Int(ceil(steps/2))
    if check
        prob =  min(1.0, steps*graph_data(bud).budbreak)
        return rand() < prob
    else
        error("No meristem found in branch")
    end
end

# ╔═╡ fb31e01e-efda-4e72-828c-efd5a286d673
branch_rule = Rule(
	Bud,
	lhs = prob_break,
	rhs = bud -> BudNode() + Internode() + Meristem()
)

# ╔═╡ 421e7e84-3c51-4e64-adc4-2a98a4536844
axiom = Internode() + Meristem();

# ╔═╡ 04d3fc8a-2ab7-4085-8fd6-4ab3f09eeebe
tree = Graph(axiom = axiom, rules = (meristem_rule, branch_rule), data = treeparams());

# ╔═╡ 1da39117-7f6d-4907-aad7-bcc0d0fc4be4
getInternode = Query(Internode);

# ╔═╡ 1bb9f437-44fb-4a7c-9eee-2730d1487e87
function elongate!(tree, query)
    for x in apply(tree, query)
        x.D = x.D .* [1.0 , 1.0 + data(tree).growth]
    end
end

# ╔═╡ 0fda99f0-ac5e-4080-bc44-3378ef41c856
function growth!(tree, query)
    elongate!(tree, query)
    rewrite!(tree)
end

# ╔═╡ 7e7fdcd0-e59e-471a-817a-87a8403440e4
function simulate(tree, query, nsteps)
    new_tree = deepcopy(tree)
    for i in 1:nsteps
        growth!(new_tree, query)
    end
    return new_tree
end

# ╔═╡ e6bdd4a1-135d-4c1b-be8e-24d5e1a3f3b1
newtree = simulate(tree, getInternode, 6);

# ╔═╡ e0ca99d2-84e1-47e3-93bc-324465f5abd3
render(Mesh(newtree))

# ╔═╡ 8c376ccb-1a52-4698-8b0f-7cb77911e8f0
md"### Ray tracing"

# ╔═╡ 23e2b273-3a58-49d3-b12e-1f96c6a497cd
md"We create the light sources corresponding to a realistic sky, using the functionality from `SkyDomes.jl`. The `VirtualPlantLab` tutorial defines an average sky for an entire day, which is possible because they model carbon dynamics on the day scale. We deviate here and create a sky corresponding with a single point in time, as we simulate plant growth continuously throughout the day."

# ╔═╡ b83835d2-9ace-4275-90b8-d1e70be8f598
waveband_conversion(Itype = :diffuse, waveband = :PAR, mode = :power)

# ╔═╡ 8e7fa801-719e-4860-abfe-39b7e4ff47a6
function create_sky(day_fraction; mesh, lat = 52.0*π/180.0, DOY = 182)
    # Compute solar irradiance
    sky_light = clear_sky(lat = lat, DOY = DOY, f = day_fraction) # W/m2
    # Conversion factors to PAR for direct and diffuse irradiance
    f_dir = waveband_conversion(Itype = :direct,  waveband = :PAR, mode = :power)
    f_dif = waveband_conversion(Itype = :diffuse, waveband = :PAR, mode = :power)
    # Actual irradiance per waveband
    Idir_PAR = f_dir * sky_light[:Idir]
    Idif_PAR = f_dif * sky_light[:Idif]
    # Create the dome of diffuse light
    dome = sky(
		mesh,
		Idir = 0.0, ## No direct solar radiation
		Idif = Idif_PAR, ## Daily Diffuse solar radiation
		nrays_dif = 1_000_000, ## Total number of rays for diffuse solar radiation
		sky_model = StandardSky, ## Angular distribution of solar radiation
		dome_method = equal_solid_angles, # Discretization of the sky dome
		ntheta = 9, ## Number of discretization steps in the zenith angle
		nphi = 12 ## Number of discretization steps in the azimuth angle
	) 
	# Add direct source
	append!(dome, sky(
			mesh, Idir = Idir_PAR, nrays_dir = 1_000_000,
			Idif = 0.0, nrays_diff = 0, theta_dir = sky_light[:theta], phi_dir = sky_light[:phi]
		)
	)
	
    return dome
end

# ╔═╡ c57ecaee-62bb-4ce0-b7e8-cfad1328bdb0
md"Create a geometry for the soil surface for light reflection."

# ╔═╡ bc88b630-43d2-4a02-a709-a80799abb837
function create_soil()
    soil = Rectangle(length = 21.0, width = 21.0)
    rotatey!(soil, π/2) ## To put it in the XY plane
	
    return soil
end

# ╔═╡ 3194f692-3150-40aa-99cd-769ca95e50a4
md"Create a geometry for the entire scene, consisting of our tree and the soil surface. We also add a material for the soil surface."

# ╔═╡ eaee419f-e008-40c4-aad9-7dfd5d6fb01a
function create_scene(tree)
    mesh = Mesh(tree)
    soil = create_soil()
    soil_material = Lambertian(τ = 0.0, ρ = 0.21)
    add!(mesh, soil, materials = soil_material)
    
    return mesh
end

# ╔═╡ d20064f5-9d30-48eb-952f-8a5c667cd8e3
md"Run the ray tracer. This consists of creating the geometry of our scene, creating our light sources, defining the ray tracer with settings of choice and finally running it with `trace!`."

# ╔═╡ 4a4378a8-b15e-4402-80cf-24b2508d2375
function run_raytracer!(tree; day_fraction = 0.5)
	mesh = create_scene(tree)
	# the directional light sources created with `SkyDomes.jl` require an "accelerated" mesh to run
	accmesh = accelerate(mesh, acceleration = BVH, rule = SAH{3}(5, 10)) 
	sources = create_sky(day_fraction; mesh = accmesh)
	# note we set `nx` and `ny` to 0: these are grid cloning parameters used to minimize the boundary effects of only simulating a part of a forest. however, as we simulate only a single tree, this is not applicable here.
	settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 0, ny = 0, parallel = true)
	raytracer = RayTracer(mesh, sources; settings)
	trace!(raytracer)
	
	return nothing
end

# ╔═╡ 448b6519-a7a7-4b01-926f-573a0153d42b
md"Running the ray tracer as we solve our differential equations is not computationally feasible. Instead, we run the ray tracer beforehand for a number of times throughout the day and define an interpolation function for every leaf that calculates the incoming PAR for a given time based on these PAR samples."

# ╔═╡ bdcb4eb2-d2ea-4edf-bbde-fa65cb4528a6
# ╠═╡ show_logs = false
function precalculate_PAR!(tree; Δf = 0.05)
	# remove potential existing PAR samples
	for node in getnodes(tree)
		if getstructmod(node) == :Leaf
			empty!(data(node).PAR_samples)
		end
	end

	# generate new PAR samples and divide by surface area
	for day_fraction in 0.0:Δf:1.0
		run_raytracer!(tree; day_fraction)
		for node in getnodes(tree)
			if getstructmod(node) == :Leaf
				PAR = data(node).mat |> power |> only
				PAR_flux = PAR / 
					(surface_area(Cuboid(), data(node).D) * 1e-4) # cm^2 to m^2
				push!(data(node).PAR_samples, PAR_flux)
			end
		end
	end

	# interpolate between samples
	for node in getnodes(tree)
		if getstructmod(node) == :Leaf
			PAR_interpolation = LinearInterpolation(
				data(node).PAR_samples,
				0.0:Δf:1.0
			);
			
			PAR_func(t, t_sunrise, t_sunset) = 
				(t % 24 - t_sunrise) / (t_sunset - t_sunrise) |>
				f -> clamp(f, 0, 1) |>
				PAR_interpolation
			
			data(node).PAR_func = PAR_func
		end
	end
	
	return nothing
end

# ╔═╡ d03ea181-c1ed-436c-af2c-4546d58287e5
md"#### Running the ray tracer"

# ╔═╡ 57d62420-46a3-4eb9-b49b-7ab16b48f4a0
md"Let's run the ray tracer and visualise the incoming PAR over time for all leaves."

# ╔═╡ 6848349a-51c4-4f89-98df-16784ea140b6
# ╠═╡ show_logs = false
precalculate_PAR!(newtree)

# ╔═╡ 4962391e-f21d-494b-a3ed-c2f938f25c2f
md"""
We can see differentiate between three types of leaf based on the plot below:
- Leaves that don't get hit by direct sunlight and only receive a small amount of PAR from diffuse sunlight.
- Leaves that get a lot of sunlight around noon, corresponding to leaves in direct sunlight that are angled mostly up.
- Leaves that get a lot of sunlight in the morning and evening, corresponding to leaves in direct sunlight that are angle mostly sideward.
"""

# ╔═╡ b2fcfbb3-968c-4ace-b344-484134c72df4
begin
	plot(xlims = (0.0, 24.0), xlabel = "Time (h)", ylabel = "PAR (W/m²)")
	for node in getnodes(newtree)
		if getstructmod(node) == :Leaf
			plot!(t -> data(node).PAR_func(t, 8.0, 20.0))
		end
	end
	plot!()
end

# ╔═╡ 7edd2626-a111-45bb-a62c-a526d8698f86
md"### Defining the environment"

# ╔═╡ 868e3bb7-2c0d-41a5-83ce-82e2e955e45a
md"""
The environment is not considered in the `VirtualPlantLab` tutorials, but we will define it here as we require it for the simulation of water dynamics.
"""

# ╔═╡ 27305ba4-1921-46b6-b088-b8e8f061e55a
Base.@kwdef struct rootparams
	split_prob::Float64 = 0.1
	Δα_μ::Float64 = pi / 10
	Δα_σ::Float64 = 0.1
end

# ╔═╡ 3ff77fe1-7d7a-46e7-9875-633927826472
struct Root <: VirtualPlantLab.Node
	D::Vector{Float64}
	α::Float64
	coords::Vector{Float64}
end

# ╔═╡ 6598fa4f-cee1-4476-b1f7-bb5e05efceea
function move(root, Δα_μ)
	Δcoords = [cos(data(root).α), sin(data(root).α)] * data(root).D[2]
	moved_root = Root(
		data(root).D,
		data(root).α + Δα_μ + graph_data(root).Δα_σ * randn(),
		data(root).coords + Δcoords
	)
	
	return moved_root
end

# ╔═╡ 809baa85-b47c-49c7-af79-1d9e126706da
root_rule = Rule(
	Root,
	lhs = root -> !has_children(root),
	rhs = root -> (
		rand() < graph_data(root).split_prob ? 
		Root(data(root).D, data(root).α, data(root).coords) + 
			(
				move(root, -graph_data(root).Δα_μ),
				move(root, graph_data(root).Δα_μ)
			) : 
		Root(data(root).D, data(root).α, data(root).coords) + 
			move(root, 0.0)
	)
)

# ╔═╡ 829d6e13-39da-4abc-ba42-eed11f9c9f70
begin
	root = Graph(axiom = Root([0.5, 10.0], -pi/2, [0.0, 0.0]), rules = (root_rule), data = rootparams())
	for _ in 1:20
		rewrite!(root)
	end
end;

# ╔═╡ 23c352ec-2a37-429e-ab6c-faba029398c0
plotstructure(root)

# ╔═╡ 5e5c75de-67a4-4578-b5c8-2537874c7043
[data(node).coords for node in getnodes(root)] |>
	x -> Plots.scatter(first.(x), last.(x))

# ╔═╡ b420e466-81e7-4328-8e1d-545133d468d9
voxel_size = 20

# ╔═╡ 99bb4507-1346-4cf3-b20d-0183a30d209d
plotstructure(root_graph)

# ╔═╡ 1ced3c28-55be-4f4e-99d0-c38f6cc79396
struct Soil <: VirtualPlantLab.Node 
	id
	x
	z
end

# ╔═╡ 8a1168a1-a678-4a78-a0b7-e1817ef58152
soil_graph = [Soil(x+5*z, x, z) for x in -2:2, z in 0:-1:-4]

# ╔═╡ b7bdc909-bfb5-4618-9702-74c77f019b2c
plotstructure(soil_graph)

# ╔═╡ 082576ba-0932-432b-b173-c844fb57bfc0
struct Air <: VirtualPlantLab.Node end

# ╔═╡ 3b249efd-55cb-4a65-a245-07296809c6b6
graphs = [newtree, soil_graph, Air()];

# ╔═╡ 6461c5f0-5860-433d-8627-53a14dd342e4
intergraph_connections = [(1, 2) => (getnodes(newtree)[1], soil_graph[1]), (1, 3) => (:Leaf, :Air)];

# ╔═╡ 52aa8d2e-16f5-49fc-86fd-0f1c3a325b3d
plantstructure = PlantStructure(graphs, intergraph_connections);

# ╔═╡ ca063cc1-cd0e-48de-b65b-0e42c0f0df65
plotstructure(plantstructure)

# ╔═╡ ac0e93c0-2662-4ea4-bcc3-3067aec41d28
md"## Functional definition"

# ╔═╡ e4c15fd1-c4f2-4cfc-a4de-8f6af50ec0de
md"""
For the definition of the functional processes, we return to `PlantModules.jl` territory. Our structural modules are assigned the same classic sets of functional modules as in previous tutorials, with two exceptions:
- We don't simulate hydraulics for the meristem, which is only a thin layer of cells.
- We simulate light interception and carbon assimilation in the leaves, with a guest appearance of the [`PlantBioPhysics.jl`](https://vezy.github.io/PlantBiophysics.jl/stable/) package (also from the `VirtualPlantLab` ecosystem) for the calculation of carbon assimilation. 
"""

# ╔═╡ 338fcbe4-f8ed-4a5d-a496-34e8199c35fb
md"### Functional modules"

# ╔═╡ 767cb212-1e7d-4909-8be2-d873962c6294
import PlantModules: t, d

# ╔═╡ 8e804ec6-fd1b-451b-b08a-d222a982b4b1
md"""
Our meristem is part of the graph, so it needs to respect the `PlantModules.multi_connection_eqs` function, which defines how nodes interact with their neighbours. The default value for this function simply states that the net water inflow $ΣF$ of a node equals the sum of the water flows $F$ between the node and each of its neighbors. Therefore, if we want the meristem to be hydraulically inactive, we need to set both its water inflow $ΣF$ and each of its water flows $F$ to $0$.
"""

# ╔═╡ 9232633a-5df0-4825-b6a9-d7efd5e036d7
function inactive_hydraulic_module(; name)
	@variables ΣF(t)
	eqs = [d(ΣF) ~ 0]

	return System(eqs, t; name)
end

# ╔═╡ 396c9109-d930-4f89-b280-ab23b918391f
function inactive_hydraulic_connection(; name)
    @variables F(t) = 0.0

    eqs = [d(F) ~ 0]
	# `get_connection_eqset` defines how variables of the connection relate to variables of two nodes it connects, but as we use no variables here that correspond with those from the nodes, the set is simply empty 
    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = []
	
    return System(eqs, t; name), get_connection_eqset
end

# ╔═╡ 20bcd5f6-fc5d-4ef8-9ea6-ebcaf2ae2bb9
md"""
Based on our ray tracer, we can calculate the incoming PAR for every leaf at any time of day. However, we still need to model how incoming PAR relates to carbon assimilation. We will use the established carbon assimilation and stomatal conducatance models available in `PlantBiophysics.jl` for this.
"""

# ╔═╡ 67223c01-0a1a-43cf-995b-429f72353536
function get_assimilation_rate(PAR_flux, T)
	Kelvin_to_C = -273.15
	meteo = Atmosphere(
		T = T + Kelvin_to_C, Wind = 1.0, P = 101.3, Rh = 0.65
	)
	m = ModelList(
		Fvcb(), # calculate CO2 assimilation rate
		Medlyn(0.03, 0.92), # calculate stomatal conductance, see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2486.2010.02375.x
		status = (Tₗ = meteo[:T], Cₛ = meteo[:Cₐ], Dₗ = meteo[:VPD], RI_PAR_f = meteo[:Ri_PAR_f], aPPFD = PAR_flux)
	)
	run!(m, meteo)
    assimilation_rate = m[:A][1] # extract result of the first (and only) timestep
	assimilation_rate = max(0, assimilation_rate) # set negative values to 0
	return assimilation_rate
end

# ╔═╡ d577f21d-bc94-4fc5-a693-6a89646eb93a
md"The above function is too complex for `ModelingToolkit.jl` to perform symbolic transformations on, and will throw an error if we use it as-is. We can exclude it from being considered for symbolic transformations using `@register_symbolic`"

# ╔═╡ 87be6a5e-ca16-4aea-aef4-8911dbb604d4
@register_symbolic get_assimilation_rate(PAR_flux, T)

# ╔═╡ 7982a512-5f3a-4c10-bdd9-6c3d6063d544
md"""
Now we can define the photosynthesis module itself. It is rather straightforward, as we have already defined our functions for getting incoming PAR based on the time of day and carbon assimilation based on the incoming PAR. Note that to access the `PAR_func` defined in the graph nodes we treat `PAR_func` as any other parameter: we include it in the inputs of our functional module and assign it a default parameter value later on. This way, our functional module will grab the value for `PAR_func` defined in the graph nodes during model creation.
"""

# ╔═╡ 4b89326d-1896-4aa5-b46d-65f2e1281255
function photosynthesis_module(; name, T, M, shape, M_c, PAR_func, t_sunrise, t_sunset)
	@constants (
		uc = (10^-6 * 10^-4 * 60^2), [description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)"],
		    # the output from PlantBiophysics.jl is in different units than we use for our ODEs, so we need to change this
	)
	@parameters (
		T = T, [description = "Temperature"], # K
        M_c = M_c, [description = "Rate of carbon consumption"], # hr^-1
	)
	@variables (
        M(t) = M, [description = "Osmotically active metabolite content"], 
			# mol / cm^3
		PF(t), [description = "Incoming PAR flux"], # W / m^2
		A(t), [description = "Carbon assimilation rate"], # mol / cm^2 / hr
		A_V(t), [description = "Volumetric carbon assimilation rate"], # mol / cm^3 / hr
		D(t)[1:getdimensionality(shape)], [description = "Dimensions of compartment"], # cm
    )

    eqs = [
		PF ~ PAR_func(t, t_sunrise, t_sunset)
		A ~ uc * get_assimilation_rate(PF, T)
		A_V ~ A * surface_area(shape, D) / PlantModules.volume(shape, D)
        d(M) ~ A_V - M_c*M
    ]
    return System(eqs, t; name, checks = false)
end

# ╔═╡ 62257b82-ac54-4360-951c-9192edb619a3
md"### Coupling"

# ╔═╡ ee971038-0d51-4fcd-9b18-c580d727265a
md"""
We couple structural and functional modules as discussed at the start of this section.
"""

# ╔═╡ 251fd6eb-627b-4560-b692-d02bef8f089c
module_coupling = Dict(
	:Meristem => [inactive_hydraulic_module],
	:Bud => [hydraulic_module, constant_carbon_module, K_module],
    :Node => [hydraulic_module, constant_carbon_module, K_module],
	:BudNode => [hydraulic_module, constant_carbon_module, K_module],
	:Internode => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, photosynthesis_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
);

# ╔═╡ 812ec2a4-b3ec-4240-a7bd-b80ddea7a748
connecting_modules = Dict(
	(:Soil, :Internode) => constant_hydraulic_connection,
	(:Internode, :Node) => hydraulic_connection,
	(:Node, :Bud) => hydraulic_connection,
	(:Node, :BudNode) => hydraulic_connection,
	(:BudNode, :Internode) => hydraulic_connection,
	(:Node, :Leaf) => hydraulic_connection,
	(:Internode, :Meristem) => inactive_hydraulic_connection,
	(:Leaf, :Air) => daynight_hydraulic_connection
);

# ╔═╡ 1f90f558-272b-41b6-a4ba-c7bbcea3cf92
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ 8f986516-75a7-46b2-86f5-4836d8b54193
md"### Parameters"

# ╔═╡ 3ef6e16e-10bc-45f4-878e-eef97d0d9829
md"""
Parameter specification follows the usual steps of assigning the correct shapes to non-cylindrical structural modules and setting the water capacity `W_max` of our soil, initial relative water content `W_r` of the air, and hydraulic conductivity `K` of the air. Additionally, we also change the value for our rate of carbon consumption `M_c`, used in our photosynthesis module. Carbon consumption is modelled as simple exponential decay. This is an extremely simplistic model, but as realistic carbon consumption is not the goal of this tutorial, it works well enough here if we adjust the rate a little.
"""

# ╔═╡ 6452a19d-fe5a-48b6-8ec8-d8759407a4dc
default_changes = Dict([:PAR_func => zero])

# ╔═╡ 1d8aba86-aba8-4567-a3f6-648d3b74e8e3
module_defaults = Dict(
	:Bud => Dict(:shape => PlantModules.Sphere()),
	:Leaf => Dict(:shape => PlantModules.Cuboid(), :M_c => 0.5),
	:Soil => Dict(:W_max => 1e4),
	:Air => Dict(:W_r => 0.6, :K => 1e-3)
);

# ╔═╡ 2f54285c-50cf-49e4-b2b0-cf9aa5fdf585
plantparams = PlantParameters(; default_changes, module_defaults);

# ╔═╡ 25c85fd1-1921-4707-8185-735751c0d914
md"## Creating the system"

# ╔═╡ 10236da6-fcef-49bf-9c3c-21f75e83f554
md"""
Finally, we generate and run the system. We have reached pure `PlantModules.jl`/`ModelingToolkit.jl`/`DifferentialEquations.jl` nirvana and therefore all steps are the same as always.
"""

# ╔═╡ 8f635abb-c9ea-45ca-81ae-e679e0ba923d
# ╠═╡ disabled = true
#=╠═╡
system = generate_system(plantstructure, plantcoupling, plantparams);
  ╠═╡ =#

# ╔═╡ 9b90ddc4-cee5-4063-9f83-ecbff8c7596c
tspan = (0.0, 2*24.0)

# ╔═╡ 7b8bbd0f-c650-49ef-b8d9-c42cd6d8da9e
#=╠═╡
prob = ODEProblem(system, [], tspan, sparse = true);
  ╠═╡ =#

# ╔═╡ 8363c723-74e0-49bf-88ae-e164919f4681
#=╠═╡
sol = solve(prob, FBDF());
  ╠═╡ =#

# ╔═╡ 10172f86-dc0b-45df-af00-a27931dca564
#=╠═╡
plotgraph(sol, plantstructure, varname = :PF, structmod = :Leaf)
  ╠═╡ =#

# ╔═╡ 4110fa63-ba4c-4aab-b9ba-ac6a4e4bd985
#=╠═╡
plotgraph(sol, plantstructure, varname = :A_V, structmod = :Leaf)
  ╠═╡ =#

# ╔═╡ d193b1aa-5498-4a4c-b2da-2a7bb2e72268
#=╠═╡
plotgraph(sol, plantstructure, varname = :M, structmod = :Leaf)
  ╠═╡ =#

# ╔═╡ 64724bb9-05c2-448c-8f72-edb88edbce45
#=╠═╡
plotgraph(sol, plantstructure, varname = :V, structmod = :Internode)
  ╠═╡ =#

# ╔═╡ 8743b00d-5699-4201-9a23-0061552495a0
#=╠═╡
plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Internode, :Bud, :Soil])
  ╠═╡ =#

# ╔═╡ 6827996d-d0c6-4d00-a13a-32495b2bf6c8
#=╠═╡
plotgraph(sol, plantstructure, varname = :P, structmod = [:Internode, :Leaf])
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─2dfdb97f-361c-47bd-b65a-41b02bc8bc57
# ╟─87a39d4c-6f22-447e-bd45-1f133b5bb7b3
# ╟─e52a8761-f703-4dec-b02b-62d3cc831b4f
# ╠═16e51c70-fe21-40c1-98f0-404254a71b1f
# ╠═ee4fa627-a613-4128-9fb5-5b336a78f18d
# ╠═6dae54ca-aa7f-4504-b985-195099162109
# ╠═813b229f-c17e-4a10-9945-dc9ad9066724
# ╠═b1f03d85-2fa4-4746-a5d8-606671d8375e
# ╠═2ea18a37-041a-48a4-9af0-d3e2bc109c63
# ╠═dcb943cf-762f-41dd-9231-e07e9d9cda7e
# ╠═5f494ba3-8c68-4728-9627-3dac9fe7fcd6
# ╠═79690228-9815-4bf4-babf-329050d9aca0
# ╠═d40434c5-4fda-4f2b-a37b-28137924a236
# ╠═986648b3-6c8a-465a-87a0-e0ce5fbbefa8
# ╠═6b1b3872-8dec-45fa-99d9-9a6903e2170a
# ╟─a36a2a5b-de0b-46a5-a97c-d0c6b1f2f9f4
# ╟─25160d35-aa6b-4ed5-9199-a54d882dbfc9
# ╟─d3fdc6b7-e69d-425e-ae3b-99e9497b8cb5
# ╟─54eb2298-38cf-49dc-8ca9-44872f0d2375
# ╟─83fadc9e-74a5-4ddd-915a-0b8692c9280f
# ╟─3632551c-2318-4bdf-b813-436c5da7dc68
# ╠═b27f3786-84b0-469f-ab0c-a60f02030efd
# ╠═b99eb3ee-26dc-4396-9a61-aba45f7462cf
# ╠═43bb6044-e567-4f20-8be0-5108f5deead8
# ╠═4ea2dfcf-cc13-46f1-8f70-f42af7a43a09
# ╠═fc109086-c42c-474d-8ab0-4add40a07eed
# ╠═2ff96078-4a91-427b-8c39-689d47600c05
# ╟─5a1e5b8a-be42-45e8-9abc-5cecd9dc83ef
# ╟─da673a6c-c661-4d00-8e42-365250cb00d3
# ╠═1f2b0dc8-84db-4612-841a-6ddaf08345e5
# ╠═e062618b-3913-4dd8-8ef4-fdcba785fa61
# ╠═d6fa9897-05ae-4388-a299-c08cb62c1414
# ╟─2c04e08b-1006-4e88-b754-1297eed4dd9c
# ╟─f69abf27-2c51-475b-8344-97361ace66d6
# ╠═120217e8-c708-4a49-8f2c-1191690ce0db
# ╠═eb87190d-26e5-4e37-bb09-15b1e89da868
# ╠═d0e32e1a-733e-4791-82ec-c4a20745a351
# ╠═fb31e01e-efda-4e72-828c-efd5a286d673
# ╠═421e7e84-3c51-4e64-adc4-2a98a4536844
# ╠═04d3fc8a-2ab7-4085-8fd6-4ab3f09eeebe
# ╠═1da39117-7f6d-4907-aad7-bcc0d0fc4be4
# ╠═1bb9f437-44fb-4a7c-9eee-2730d1487e87
# ╠═0fda99f0-ac5e-4080-bc44-3378ef41c856
# ╠═7e7fdcd0-e59e-471a-817a-87a8403440e4
# ╠═e6bdd4a1-135d-4c1b-be8e-24d5e1a3f3b1
# ╠═e0ca99d2-84e1-47e3-93bc-324465f5abd3
# ╟─8c376ccb-1a52-4698-8b0f-7cb77911e8f0
# ╟─23e2b273-3a58-49d3-b12e-1f96c6a497cd
# ╠═b83835d2-9ace-4275-90b8-d1e70be8f598
# ╠═8e7fa801-719e-4860-abfe-39b7e4ff47a6
# ╟─c57ecaee-62bb-4ce0-b7e8-cfad1328bdb0
# ╠═bc88b630-43d2-4a02-a709-a80799abb837
# ╟─3194f692-3150-40aa-99cd-769ca95e50a4
# ╠═eaee419f-e008-40c4-aad9-7dfd5d6fb01a
# ╟─d20064f5-9d30-48eb-952f-8a5c667cd8e3
# ╠═4a4378a8-b15e-4402-80cf-24b2508d2375
# ╟─448b6519-a7a7-4b01-926f-573a0153d42b
# ╠═bdcb4eb2-d2ea-4edf-bbde-fa65cb4528a6
# ╟─d03ea181-c1ed-436c-af2c-4546d58287e5
# ╟─57d62420-46a3-4eb9-b49b-7ab16b48f4a0
# ╠═6848349a-51c4-4f89-98df-16784ea140b6
# ╟─4962391e-f21d-494b-a3ed-c2f938f25c2f
# ╠═b2fcfbb3-968c-4ace-b344-484134c72df4
# ╟─7edd2626-a111-45bb-a62c-a526d8698f86
# ╟─868e3bb7-2c0d-41a5-83ce-82e2e955e45a
# ╠═27305ba4-1921-46b6-b088-b8e8f061e55a
# ╠═3ff77fe1-7d7a-46e7-9875-633927826472
# ╠═6598fa4f-cee1-4476-b1f7-bb5e05efceea
# ╠═809baa85-b47c-49c7-af79-1d9e126706da
# ╠═829d6e13-39da-4abc-ba42-eed11f9c9f70
# ╠═23c352ec-2a37-429e-ab6c-faba029398c0
# ╠═5e5c75de-67a4-4578-b5c8-2537874c7043
# ╠═b420e466-81e7-4328-8e1d-545133d468d9
# ╠═99bb4507-1346-4cf3-b20d-0183a30d209d
# ╠═1ced3c28-55be-4f4e-99d0-c38f6cc79396
# ╠═8a1168a1-a678-4a78-a0b7-e1817ef58152
# ╠═b7bdc909-bfb5-4618-9702-74c77f019b2c
# ╠═082576ba-0932-432b-b173-c844fb57bfc0
# ╠═3b249efd-55cb-4a65-a245-07296809c6b6
# ╠═6461c5f0-5860-433d-8627-53a14dd342e4
# ╠═52aa8d2e-16f5-49fc-86fd-0f1c3a325b3d
# ╠═ca063cc1-cd0e-48de-b65b-0e42c0f0df65
# ╟─ac0e93c0-2662-4ea4-bcc3-3067aec41d28
# ╟─e4c15fd1-c4f2-4cfc-a4de-8f6af50ec0de
# ╟─338fcbe4-f8ed-4a5d-a496-34e8199c35fb
# ╠═767cb212-1e7d-4909-8be2-d873962c6294
# ╟─8e804ec6-fd1b-451b-b08a-d222a982b4b1
# ╠═9232633a-5df0-4825-b6a9-d7efd5e036d7
# ╠═396c9109-d930-4f89-b280-ab23b918391f
# ╟─20bcd5f6-fc5d-4ef8-9ea6-ebcaf2ae2bb9
# ╠═67223c01-0a1a-43cf-995b-429f72353536
# ╟─d577f21d-bc94-4fc5-a693-6a89646eb93a
# ╠═87be6a5e-ca16-4aea-aef4-8911dbb604d4
# ╟─7982a512-5f3a-4c10-bdd9-6c3d6063d544
# ╠═4b89326d-1896-4aa5-b46d-65f2e1281255
# ╟─62257b82-ac54-4360-951c-9192edb619a3
# ╟─ee971038-0d51-4fcd-9b18-c580d727265a
# ╠═251fd6eb-627b-4560-b692-d02bef8f089c
# ╠═812ec2a4-b3ec-4240-a7bd-b80ddea7a748
# ╠═1f90f558-272b-41b6-a4ba-c7bbcea3cf92
# ╟─8f986516-75a7-46b2-86f5-4836d8b54193
# ╟─3ef6e16e-10bc-45f4-878e-eef97d0d9829
# ╠═6452a19d-fe5a-48b6-8ec8-d8759407a4dc
# ╠═1d8aba86-aba8-4567-a3f6-648d3b74e8e3
# ╠═2f54285c-50cf-49e4-b2b0-cf9aa5fdf585
# ╟─25c85fd1-1921-4707-8185-735751c0d914
# ╟─10236da6-fcef-49bf-9c3c-21f75e83f554
# ╠═8f635abb-c9ea-45ca-81ae-e679e0ba923d
# ╠═9b90ddc4-cee5-4063-9f83-ecbff8c7596c
# ╠═7b8bbd0f-c650-49ef-b8d9-c42cd6d8da9e
# ╠═8363c723-74e0-49bf-88ae-e164919f4681
# ╠═10172f86-dc0b-45df-af00-a27931dca564
# ╠═4110fa63-ba4c-4aab-b9ba-ac6a4e4bd985
# ╠═d193b1aa-5498-4a4c-b2da-2a7bb2e72268
# ╠═64724bb9-05c2-448c-8f72-edb88edbce45
# ╠═8743b00d-5699-4201-9a23-0061552495a0
# ╠═6827996d-d0c6-4d00-a13a-32495b2bf6c8
