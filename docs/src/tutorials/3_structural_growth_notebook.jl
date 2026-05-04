### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 16e51c70-fe21-40c1-98f0-404254a71b1f
using Pkg; Pkg.activate("../..")

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
In the previous two tutorials, we simulated water dynamics for a plant given a certain, static structure. This tutorial, we will also consider the **structural** growth of a plant. As `PlantModules.jl` has no inherent structural modelling capabilities, we will instead how to integrate with the graph rewriting functionality provided by [`VirtualPlantLab.jl`](https://virtualplantlab.com/stable) to achieve this. Specifically, we will couple our package's functionality to the rewriting process by making the rewrite rules depend on the functional status of the tree.

In addition, this tutorial considers two more geometry-related functionalities typically expected of FSPMs:
- We simulate carbon dynamics in the leaves based on a ray tracer, a carbon assimilation model and a stomatal conductance model, all of which are available in `VirtualPlantLab.jl`.
- We discretize the soil into multiple compartments to simulate water transport more realistically.
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

For didactical purposes, we will base this tutorial directly on the series of tutorials mentioned above. Note that in `VirtualPlantLab.jl`'s ray tracing tutorial, light interception is used to simulate carbon dynamics based on a discrete sink-source implentation. We will deviate here and instead use `PlantModules.jl`'s functionality to model carbon dynamics coupled to water relations based on differential equations.

The tree structure used in the ray tracing tutorial has a lot of additional information for the discrete sink-source simulation of carbon dynamics. Considering we don't require this information, we will use the simpler tree structure from the [first tutorial](https://virtualplantlab.com/stable/tutorials/from_tree_forest/tree/) in the series.
"""

# ╔═╡ 33038274-3206-4433-9f6d-79dccff654cb
md"### Shoots"

# ╔═╡ 83fadc9e-74a5-4ddd-915a-0b8692c9280f
md"#### Structural modules"

# ╔═╡ 3632551c-2318-4bdf-b813-436c5da7dc68
md"""
The tree consists of the classic set of internodes, nodes and leaves, as well as buds that grow into branches, corresponding nodes left behind by the buds, and meristems that grow into additional phytomers. 

Some things of note for this tutorial:
- All plant parts that are plotted need to have their dimensions defined;
- All plant parts that intercept light in the ray tracing additionally need to have a material defined (see the VPL docs);
- Plant parts that are part of a rewriting rule depending on the plant's functional status, should have those functional variables defined and be defined to be mutable.
"""

# ╔═╡ b27f3786-84b0-469f-ab0c-a60f02030efd
Base.@kwdef struct Meristem <: VirtualPlantLab.Node end

# ╔═╡ b99eb3ee-26dc-4396-9a61-aba45f7462cf
Base.@kwdef mutable struct Bud <: VirtualPlantLab.Node
	D::Vector{Float64} = [0.5]
	W::Float64 = PlantModules.volume(PlantModules.Sphere(), [0.5])
end

# ╔═╡ 43bb6044-e567-4f20-8be0-5108f5deead8
Base.@kwdef struct Node <: VirtualPlantLab.Node end

# ╔═╡ 4ea2dfcf-cc13-46f1-8f70-f42af7a43a09
Base.@kwdef struct BudNode <: VirtualPlantLab.Node end

# ╔═╡ fc109086-c42c-474d-8ab0-4add40a07eed
Base.@kwdef mutable struct Internode <: VirtualPlantLab.Node
	D::Vector{Float64} = [0.5, 3.0] # newly grown internodes start small
	mat::Lambertian{1} = Lambertian(τ = 0.05, ρ = 0.1) # internodes require a material as they interact with light
end

# ╔═╡ 2ff96078-4a91-427b-8c39-689d47600c05
Base.@kwdef mutable struct Leaf <: VirtualPlantLab.Node
	D::Vector{Float64} = [5, 3, 0.05]
	mat::Lambertian{1} = Lambertian(τ = 0.05, ρ = 0.1) # leaves also interact with light
	PAR_samples::Vector{Float64} = Float64[]
	PAR_func::Function = (x...) -> error("A leaf node exists with no defined PAR function.")
end

# ╔═╡ 5a1e5b8a-be42-45e8-9abc-5cecd9dc83ef
md"#### Geometry functions"

# ╔═╡ da673a6c-c661-4d00-8e42-365250cb00d3
md"Define a `feed!` method for every structural module that needs a defined geometry, i.e. that needs to be plotted or interacts with light."

# ╔═╡ 1f2b0dc8-84db-4612-841a-6ddaf08345e5
function VirtualPlantLab.feed!(turtle::Turtle, i::Internode, vars)
    # Rotate turtle around the head to implement elliptical phyllotaxis
    rh!(turtle, vars.phyllotaxis)
	# remember to divide by 100 to convert cm to m
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
md"#### Structural growth rules"

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
	W_burst::Float64 = 1.5 # water content at which bud is guaranteed to burst
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
    prob = data(bud).W / graph_data(bud).W_burst
	return rand() < prob
end

# ╔═╡ fb31e01e-efda-4e72-828c-efd5a286d673
branch_rule = Rule(
	Bud,
	lhs = prob_break,
	rhs = bud -> BudNode() + Internode() + Meristem()
)

# ╔═╡ e6bdd4a1-135d-4c1b-be8e-24d5e1a3f3b1
shoot_graph = Graph(
	axiom = Internode(D = [0.5, 10.0]) + Node() + (Bud(), Leaf()) + Internode() + Meristem(),
	rules = (meristem_rule, branch_rule), data = treeparams()
);

# ╔═╡ e0ca99d2-84e1-47e3-93bc-324465f5abd3
render(Mesh(shoot_graph))

# ╔═╡ dd6debaf-62d4-43da-88a9-bebbb303d095
md"### Roots"

# ╔═╡ 432215ea-eea5-4edb-a21f-4d9909b2c4c4
md"""
We model the roots as consisting of a single structural module. At every rewriting step, it either elongates or splits in two with a set probability. As we will have to connect the root segments to the correct soil compartment later on, it's important to track the positions of the root segments. We include positional information in our model by assigning a growing direction to all root segments. When new root segments are produced, their location follows this growing direction with a certain random perturbation.
"""

# ╔═╡ 27305ba4-1921-46b6-b088-b8e8f061e55a
Base.@kwdef struct rootparams
	split_prob::Float64 = 0.2
	Δdirection_σ::Float64 = 0.2
	Δdirection_σ_branch::Float64 = 0.5
end;

# ╔═╡ 3ff77fe1-7d7a-46e7-9875-633927826472
struct Root <: VirtualPlantLab.Node
	D::Vector{Float64}
	coords::Vector{Float64}
	direction::Vector{Float64}
end

# ╔═╡ 6598fa4f-cee1-4476-b1f7-bb5e05efceea
function move(root, Δdirection_σ)
	new_D = [0.5, 3.0] # newly grown roots start small
	new_direction = data(root).direction + Δdirection_σ * randn(3) |>
		x -> x / sqrt(sum(x.^2)) # normalize to unit length
	new_coords = data(root).coords + new_D[2] * new_direction
	
	moved_root = Root(
		new_D,
		new_coords,
		new_direction
	)
	
	return moved_root
end

# ╔═╡ 809baa85-b47c-49c7-af79-1d9e126706da
root_rule = Rule(
	Root,
	lhs = root -> !has_children(root),
	rhs = root -> (
		rand() < graph_data(root).split_prob ? 
		Root(data(root).D, data(root).coords, data(root).direction) + 
			(
				move(root, graph_data(root).Δdirection_σ_branch),
				move(root, graph_data(root).Δdirection_σ_branch)
			) : 
		Root(data(root).D, data(root).coords, data(root).direction) + 
			move(root, graph_data(root).Δdirection_σ)
	)
)

# ╔═╡ 829d6e13-39da-4abc-ba42-eed11f9c9f70
root_graph = Graph(
	axiom = sum([Root([0.5, 10.0], [0.0, 0.0, i*-10.0], [0.0, 0.0, -1.0]) for i in 1:11]),
	rules = (root_rule), data = rootparams()
)

# ╔═╡ 9f4637df-e1b6-4bfc-b66c-980be1b24f19
md"""
!!! note
	Purely to demonstrate connecting with soil voxels in the next part, we start off with a long root.
"""

# ╔═╡ 23c352ec-2a37-429e-ab6c-faba029398c0
plotstructure(root_graph)

# ╔═╡ 7edd2626-a111-45bb-a62c-a526d8698f86
md"### Environment"

# ╔═╡ 868e3bb7-2c0d-41a5-83ce-82e2e955e45a
md"""
As mentioned before, we will consider a soil consisting of multiple compartments. More specifically, we will discretize the soil into voxels, or cubes, of a given size. To connect each voxel to its non-diagonal neighbours, we can simply define the graph as a 3-dimensional array of soil nodes, which will get translated into our desired graph by `PlantModules.jl`. However, as Julia arrays are no typical graph format, they store no information about the `id`, and therefore it is required to manually define the `id` inside of our structural modules. Additionally, we will include the $x$, $y$ and $z$ coordinates of the voxel's centre.
"""

# ╔═╡ 1ced3c28-55be-4f4e-99d0-c38f6cc79396
Base.@kwdef struct Soil <: VirtualPlantLab.Node 
	id::Integer
	x
	y
	z
end

# ╔═╡ e677eb3d-460f-4a42-bffb-56286c3fc6de
vs = 100.0; # voxel size

# ╔═╡ 8a1168a1-a678-4a78-a0b7-e1817ef58152
soil_graph = [
	Soil(id = x + 3*y + 9*z; x, y, z) 
	for x in vs*(-1:1), y in vs*(-1:1), z in vs*(-0.5:-1:-2.5)
];

# ╔═╡ b7bdc909-bfb5-4618-9702-74c77f019b2c
plotstructure(soil_graph)

# ╔═╡ 4263fc67-f667-4e9c-93f7-bea1ab7b9579
md"The air could be similarly divided into compartments. For simplicity, however, we will again model it as a single node."

# ╔═╡ 082576ba-0932-432b-b173-c844fb57bfc0
struct Air <: VirtualPlantLab.Node end

# ╔═╡ 8c376ccb-1a52-4698-8b0f-7cb77911e8f0
md"### Ray tracing"

# ╔═╡ 23e2b273-3a58-49d3-b12e-1f96c6a497cd
md"We create the light sources corresponding to a realistic sky, using the functionality from `SkyDomes.jl` (from the from `VirtualPlantLab.jl` ecosystem). The ray tracing tutorial we base ourselves on defines an average sky for an entire day, which is possible because they model carbon dynamics on the day scale. We deviate here and create a sky corresponding with a single point in time, as we simulate plant growth continuously throughout the day."

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
	settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 0, ny = 0,
						  parallel = true, verbose = false)
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
				0.0:Δf:1.0,
				extrapolation = ExtrapolationType.Constant # use constant boundary values for extrapolation (here always 0)
			);
			
			PAR_func(t, t_sunrise, t_sunset) = 
				(t % 24 - t_sunrise) / (t_sunset - t_sunrise) |>
				f -> clamp(f, 0.0, 1.0) |>
				PAR_interpolation
			
			data(node).PAR_func = PAR_func
		end
	end
	
	return nothing
end

# ╔═╡ d03ea181-c1ed-436c-af2c-4546d58287e5
md"#### Running the ray tracer"

# ╔═╡ 57d62420-46a3-4eb9-b49b-7ab16b48f4a0
md"Let's run the ray tracer and visualise the incoming PAR over time for our leaf."

# ╔═╡ 6848349a-51c4-4f89-98df-16784ea140b6
# ╠═╡ show_logs = false
precalculate_PAR!(shoot_graph)

# ╔═╡ b2fcfbb3-968c-4ace-b344-484134c72df4
begin
	plot(xlims = (0.0, 24.0), xlabel = "Time (h)", ylabel = "PAR (W/m²)", legend = false)
	for node in getnodes(shoot_graph)
		if getstructmod(node) == :Leaf
			plot!(t -> data(node).PAR_func(t, 8.0, 20.0))
		end
	end
	plot!()
end

# ╔═╡ ea169cfa-510f-4619-8e63-394ea0a44b09
md"### Connecting the separate graphs"

# ╔═╡ 845af412-656b-4b65-8d09-9f2bc58b3867
md"""
We finalize the structural definition by connecting all the separate graphs. The graphs are connected as follows:
- Shoots and roots: the base of the shoots is connected to the base of the root system.
- Roots and soil: each root segment is connected to the soil voxel they are physically inside of, which is done by comparing the coordinates of the root segments and soil voxels. Complicated connections like this can be defined by passing a function to `intergraph_connections` for the graphs in question that takes a node of each graph and returns whether they should be connected.
- Shoots and air: each leaf node is connected to the single air node.
- Soil and air: each node of the top soil layer is connected to the air to simulate evaporation.
"""

# ╔═╡ 3b249efd-55cb-4a65-a245-07296809c6b6
graphs = [shoot_graph, root_graph, soil_graph, Air()];

# ╔═╡ 69aa64f5-f0ed-4fa0-88bc-59141a8e42b9
function is_connected_root_soil(root, soil)
	soil_center = [soil.x, soil.y, soil.z]
	coord_inside_voxel = [
		soil_center[i]-vs/2 <= data(root).coords[i] < soil_center[i]+vs/2
		for i in eachindex(data(root).coords)
	]
	are_nodes_connected = all(coord_inside_voxel)
	return are_nodes_connected
end

# ╔═╡ 6461c5f0-5860-433d-8627-53a14dd342e4
intergraph_connections = [
	(1, 2) => (getnodes(shoot_graph)[1], getnodes(root_graph)[1]),
	(2, 3) => is_connected_root_soil,
	(1, 4) => (:Leaf, :Air),
	(3, 4) => (soil_graph[:, :, 1], :Air)
];

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

# ╔═╡ cadcf258-862b-416a-8e5e-8fb52a0e7b2a
plot(PAR -> get_assimilation_rate(PAR, 293.15), xlims = (0, 300), xlabel = "PAR (W / m²)", ylabel = "Assimilation rate (μmol / m² / s)", legend = false, title = "Assimilation rate in function of PAR")

# ╔═╡ 4b381403-6fd3-485e-a6c8-dde6dd8ccc55
md"""
When we want to use more complex functions such as this inside of our functional modules, it's often a good idea to reduce their computational cost if possible. For this example, we can simply replace the function with another linear interpolation.
"""

# ╔═╡ f4ddafdc-75e8-47ad-9250-8c8d83f58394
interpolation_range = 0:1000

# ╔═╡ 9b8b866a-b44a-4915-bce5-92cc79c70818
get_assimilation_rate_interpolation = LinearInterpolation(
	get_assimilation_rate.(interpolation_range, 293.15),
	interpolation_range,
	extrapolation = ExtrapolationType.Extension # smoothly extend interpolation for extrapolation
);

# ╔═╡ b7749f99-c0d0-489f-a0e5-05a3f29713cd
plot(par -> get_assimilation_rate_interpolation(par), xlims = (-100, 500), label = false, xlabel = "PAR (W / m²)", ylabel = "Assimilation rate (μmol / m² / s)", legend = false, title = "Inter/extrapolation of assimilation rate")

# ╔═╡ 7982a512-5f3a-4c10-bdd9-6c3d6063d544
md"""
Now we can define the photosynthesis module itself. It is rather straightforward, as we have already defined our functions for getting incoming PAR based on the time of day and carbon assimilation based on the incoming PAR. Note that to access the `PAR_func` defined in the graph nodes we treat `PAR_func` as any other parameter: we include it in the inputs of our functional module and assign it a default parameter value later on. This way, our functional module will grab the value for `PAR_func` defined in the graph nodes during model creation.
"""

# ╔═╡ 4b89326d-1896-4aa5-b46d-65f2e1281255
function photosynthesis_module(; name, M, shape, M_c, PAR_func, t_sunrise, t_sunset)
	@constants (
		uc = (10^-6 * 10^-4 * 60^2), [description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)"],
		    # the output from PlantBiophysics.jl is in different units than we use for our ODEs, so we need to change this
	)
	@parameters (
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
		A ~ uc * get_assimilation_rate_interpolation(PF)
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

	:Root => [hydraulic_module, constant_carbon_module, K_module],
	
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
);

# ╔═╡ 812ec2a4-b3ec-4240-a7bd-b80ddea7a748
connecting_modules = Dict(
	(:Soil, :Soil) => constant_hydraulic_connection,
	(:Soil, :Root) => constant_hydraulic_connection,
	(:Root, :Root) => hydraulic_connection,
	(:Root, :Internode) => hydraulic_connection,
	(:Internode, :Node) => hydraulic_connection,
	(:Node, :Bud) => hydraulic_connection,
	(:Node, :BudNode) => hydraulic_connection,
	(:BudNode, :Internode) => hydraulic_connection,
	(:Node, :Leaf) => hydraulic_connection,
	(:Internode, :Meristem) => inactive_hydraulic_connection,
	(:Leaf, :Air) => daynight_hydraulic_connection,
	(:Soil, :Air) => constant_hydraulic_connection
);

# ╔═╡ 1f90f558-272b-41b6-a4ba-c7bbcea3cf92
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ 8f986516-75a7-46b2-86f5-4836d8b54193
md"### Parameters"

# ╔═╡ 3ef6e16e-10bc-45f4-878e-eef97d0d9829
md"""
Parameter specification follows the usual steps of assigning the correct shapes and initial dimensions to the plant parts and setting the water capacity `W_max` of our soil, initial relative water content `W_r` of the air, and hydraulic conductivity `K` of the air. We also lower the hydraulic conductivity between soil and air to a more realistic value for direct evaporation from the soil.
"""

# ╔═╡ 6452a19d-fe5a-48b6-8ec8-d8759407a4dc
default_changes = Dict(
	:PAR_func => (x...) -> error("A leaf node exists with no defined PAR function.")
)

# ╔═╡ 1d8aba86-aba8-4567-a3f6-648d3b74e8e3
module_defaults = Dict(
	:Bud => Dict(:shape => PlantModules.Sphere()),
	:Leaf => Dict(:shape => PlantModules.Cuboid()),
	:Node => Dict(:D => [0.5, 0.1]),
	:BudNode => Dict(:D => [0.5, 0.1]),
	:Soil => Dict(:W_max => 1e4),
	:Air => Dict(:W_r => 0.6, :K => 1e-3)
);

# ╔═╡ fc7dbd2e-f621-4953-bd81-b565d74f0af8
connection_values = Dict(
	(:Soil, :Air) => Dict(:K => 1e-2)
)

# ╔═╡ 2f54285c-50cf-49e4-b2b0-cf9aa5fdf585
plantparams = PlantParameters(; default_changes, module_defaults, connection_values);

# ╔═╡ 25c85fd1-1921-4707-8185-735751c0d914
md"## Running the model"

# ╔═╡ 10236da6-fcef-49bf-9c3c-21f75e83f554
md"""
Finally, we generate and run the system. We will start with an initial run to generate our solution for the first day as in our previous two tutorials. For the following days, we define a function that performs the following steps:
- Update the structures of the shoot and root graphs, the former of which depends on the functional status of the plant evaluated at the end of the previous day.
- Run our ray tracer again to recalculate the incoming PAR for our updated shoot structure.
- Create an updated `PlantStructure` using the updated shoot and root graphs.
- Add the solution of the previous day's simulation to `PlantParameters`. This makes it use the final values of all variables as the initial values for the simulation of the next day.
- Generate and solve the system.
- Store our solution and plantstructure for visualisation of the results.
- Based on the day's simulation, update the water contents of bud nodes and the dimensions of the internodes and leaves in our shoot graph. These are no longer used directly by our package, as the initial values they specify are now overwritten by the final variable values specified in the previous day's solution. However, they are used in the next day's rewriting steps and ray tracing, so it is important to update them.
"""

# ╔═╡ adacf377-bbde-4367-b48c-6c4a4721dbfe
md"### Initial run"

# ╔═╡ 8f635abb-c9ea-45ca-81ae-e679e0ba923d
system = generate_system(plantstructure, plantcoupling, plantparams);

# ╔═╡ 9b90ddc4-cee5-4063-9f83-ecbff8c7596c
tspan = (0.0, 24.0);

# ╔═╡ 7b8bbd0f-c650-49ef-b8d9-c42cd6d8da9e
prob = ODEProblem(system, [], tspan, sparse = true);

# ╔═╡ 8363c723-74e0-49bf-88ae-e164919f4681
sol = solve(prob, FBDF());

# ╔═╡ 3424b05c-e792-4cb1-91ec-d27a04c8d4a2
md"### Subsequent runs"

# ╔═╡ b93903cb-24c9-4ff1-8e31-ee6cb9039c24
function run_timestep!(sols, plantstructures, shoot_graph, root_graph)
	# rewrite structures
	rewrite!(shoot_graph)
	rewrite!(root_graph)

	# run ray tracer for new leaves
	precalculate_PAR!(shoot_graph)

	# connect new structures
	graphs = [shoot_graph, root_graph, soil_graph, Air()]
	intergraph_connections = [
		(1, 2) => (getnodes(shoot_graph)[1], getnodes(root_graph)[1]),
		(2, 3) => is_connected_root_soil,
		(1, 4) => (:Leaf, :Air),
		(3, 4) => (soil_graph[:, :, 1], :Air)
	]
	plantstructure_new = PlantStructure(graphs, intergraph_connections);

	# add solution to parameters - changes initial values to previous solution's end values
	plantparams_new = PlantParameters(; default_changes, module_defaults,
									  connection_values, sol = sols[end]);

	# generate system and solve
	system_new = generate_system(plantstructure_new, plantcoupling, plantparams_new);
	prob_new = ODEProblem(system_new, [], tspan, sparse = true);
	sol_new = solve(prob_new, FBDF());

	# store data for plotting
	push!(sols, sol_new)
	push!(plantstructures, plantstructure_new)

	# update node data used in rewriting rules, plotting and ray tracing
	for node in getnodes(shoot_graph)
		if getstructmod(node) == :Bud
			var = get_subsystem_variables(
				system_new, shoot_graph, 1, :W, node
			)
			getattributes(node)[:W] = sol_new[var][end]
		elseif getstructmod(node) in [:Leaf, :Internode]
			var = get_subsystem_variables(
				system_new, shoot_graph, 1, :D, node
			)
			getattributes(node)[:D] = sol_new[var][end]
		end
	end
	
	return nothing
end

# ╔═╡ 20c98896-8ae1-4315-a517-3ecf69ca0123
days = 5

# ╔═╡ d3a39f1d-879e-4e30-9098-6cd8cc58a5af
sols = [sol];

# ╔═╡ 7aad257a-e6b6-4249-b6b2-2834e5be86d5
plantstructures = [plantstructure];

# ╔═╡ 7771dc15-9365-4722-bffc-7ba2534afaf1
for day in 2:days
	@info "Simulating day $day"
	run_timestep!(sols, plantstructures, shoot_graph, root_graph)
end

# ╔═╡ dba10e9c-1ff0-41fa-af37-86f02cd47ce2
md"Let's quickly inspect the structure of our plant at the end of the simulation period."

# ╔═╡ e33780b3-1e4c-40ac-a761-8877c79ac170
render(Mesh(shoot_graph))

# ╔═╡ 720da1fc-136d-4693-ba1c-b9e7c885e33e
plotstructure(root_graph)

# ╔═╡ c8213569-f42a-49ba-bcca-34d7d5b2b04d
plotstructure(plantstructures[end])

# ╔═╡ 7dd5d141-f6a5-467a-ac59-d766254a0e0d
md"## Results"

# ╔═╡ e27b5cce-ce85-4e21-85fe-a8757283282f
md"""
Finally, we can visualize the results of our simulation over the entire time period. We will make the following plots to show the effects of the new functionalities discussed in this tutorial
- Bud water content over time, illustrating the structural growth of our plant as new buds grow and old buds grow into branches.
- Soil water content over time, illustrating the effect of discretizing the soil into multiple compartments.
- The leaf assimilation rate over time, illustrating the different trends of incoming PAR for different leaves.

Plotting has become technically more complex now that we have multiple plantstructures and corresponding solutions, but we can simply pass them to `plotgraph` as vectors and use the plotting function as per usual. Plotting soil water content over time is more complicated, however, because we want to categorize them based on whether they are directly below the plant.
"""

# ╔═╡ 4ac9913e-5d2f-41b3-b782-4d31b441602e
plotgraph(
	sols, plantstructures, varname = :W, structmod = :Bud,
	ylabel = "Water content (g)", xlabel = "Time (h)", label = false, lw = 2,
	title = "Water content of buds", size = (800, 600), margins = 5*Plots.mm
)

# ╔═╡ d1469e76-1e73-42f5-afd3-90f1bbfb4dfa
begin
	soil_nodes = [node for node in getnodes(plantstructure) 
				  if getstructmod(node) == :Soil]
	soil_vars = get_subsystem_variables(system, plantstructure, :W, :Soil)
	soil_depths = [getattributes(node)[:z] for node in soil_nodes] |> unique |> sort
	get_label(z) = (
		z == soil_depths[1] ? "Bottom layer" :
			(z == soil_depths[2] ? "Middle layer" : "Top layer")
	)
	get_color(z) = (
		z == soil_depths[1] ? :red :
			(z == soil_depths[2] ? :orange : :blue)
	)

	center_vars = [
		soil_var 
		for (soil_var, soil_node) in zip(soil_vars, soil_nodes) 
		if getattributes(soil_node)[:x] == 0 && getattributes(soil_node)[:y] == 0
	]
	center_labels = [
		get_label(getattributes(soil_node)[:z]) 
		for soil_node in soil_nodes 
		if getattributes(soil_node)[:x] == 0 && getattributes(soil_node)[:y] == 0
	] |> x -> reshape(x, 1, :)
	center_colors = [
		get_color(getattributes(soil_node)[:z]) 
		for (soil_var, soil_node) in zip(soil_vars, soil_nodes) 
		if getattributes(soil_node)[:x] == 0 && getattributes(soil_node)[:y] == 0
	] |> x -> reshape(x, 1, :)

	border_vars = [
		soil_var
		for (soil_var, soil_node) in zip(soil_vars, soil_nodes) 
		if getattributes(soil_node)[:x] != 0 || getattributes(soil_node)[:y] != 0
	]
	border_colors = [
		get_color(getattributes(soil_node)[:z]) 
		for (soil_var, soil_node) in zip(soil_vars, soil_nodes) 
		if getattributes(soil_node)[:x] != 0 || getattributes(soil_node)[:y] != 0
	] |> x -> reshape(x, 1, :)

	xs_vec = [sol.t for sol in sols]
	for i in eachindex(xs_vec)[2:end]
		xs_vec[i] = xs_vec[i] .+ xs_vec[i-1][end]
	end
	xs = reduce(vcat, xs_vec)
	center_ys = vcat([permutedims(reduce(hcat, sol[center_vars])) for sol in sols]...)
	border_ys = vcat([permutedims(reduce(hcat, sol[border_vars])) for sol in sols]...)

	p_center = plot(
		xs, center_ys,
		label = center_labels, color = center_colors, lw = 2,
		title = "Center slice", xticks = 0:24:xs[end], ylabel = "Water content (g)"
	)
	p_border = plot(
		xs, border_ys,
		label = false, color = border_colors, lw = 2,
		title = "Border slices", xticks = 0:24:xs[end])
	
	plot(
		p_center, p_border, xlabel = "Time (h)", 
		plot_title = "Water content of soil compartments", 
		plot_titlevspan = 0.1, size = (800, 600),  margins = 5*Plots.mm
	)
end

# ╔═╡ 53ea2aff-3bcc-41ed-be7e-48b8e974f2ba
plotgraph(
	sols, plantstructures, varname = :A, structmod = :Leaf,
	ylabel = "Carbon assimilation rate (mol / cm² / h)", xlabel = "Time (h)",
	title = "Carbon assimilation rate of leaves", lw = 1.5, label = false, 
	size = (800, 600), xlims = (72, 96), xticks = 0:6:96,  margins = 5*Plots.mm
)

# ╔═╡ Cell order:
# ╟─2dfdb97f-361c-47bd-b65a-41b02bc8bc57
# ╟─87a39d4c-6f22-447e-bd45-1f133b5bb7b3
# ╟─e52a8761-f703-4dec-b02b-62d3cc831b4f
# ╠═16e51c70-fe21-40c1-98f0-404254a71b1f
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
# ╟─33038274-3206-4433-9f6d-79dccff654cb
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
# ╠═e6bdd4a1-135d-4c1b-be8e-24d5e1a3f3b1
# ╠═e0ca99d2-84e1-47e3-93bc-324465f5abd3
# ╟─dd6debaf-62d4-43da-88a9-bebbb303d095
# ╟─432215ea-eea5-4edb-a21f-4d9909b2c4c4
# ╠═27305ba4-1921-46b6-b088-b8e8f061e55a
# ╠═3ff77fe1-7d7a-46e7-9875-633927826472
# ╠═6598fa4f-cee1-4476-b1f7-bb5e05efceea
# ╠═809baa85-b47c-49c7-af79-1d9e126706da
# ╠═829d6e13-39da-4abc-ba42-eed11f9c9f70
# ╟─9f4637df-e1b6-4bfc-b66c-980be1b24f19
# ╠═23c352ec-2a37-429e-ab6c-faba029398c0
# ╟─7edd2626-a111-45bb-a62c-a526d8698f86
# ╟─868e3bb7-2c0d-41a5-83ce-82e2e955e45a
# ╠═1ced3c28-55be-4f4e-99d0-c38f6cc79396
# ╠═e677eb3d-460f-4a42-bffb-56286c3fc6de
# ╠═8a1168a1-a678-4a78-a0b7-e1817ef58152
# ╠═b7bdc909-bfb5-4618-9702-74c77f019b2c
# ╟─4263fc67-f667-4e9c-93f7-bea1ab7b9579
# ╠═082576ba-0932-432b-b173-c844fb57bfc0
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
# ╠═b2fcfbb3-968c-4ace-b344-484134c72df4
# ╟─ea169cfa-510f-4619-8e63-394ea0a44b09
# ╟─845af412-656b-4b65-8d09-9f2bc58b3867
# ╠═3b249efd-55cb-4a65-a245-07296809c6b6
# ╠═69aa64f5-f0ed-4fa0-88bc-59141a8e42b9
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
# ╟─cadcf258-862b-416a-8e5e-8fb52a0e7b2a
# ╟─4b381403-6fd3-485e-a6c8-dde6dd8ccc55
# ╠═f4ddafdc-75e8-47ad-9250-8c8d83f58394
# ╠═9b8b866a-b44a-4915-bce5-92cc79c70818
# ╟─b7749f99-c0d0-489f-a0e5-05a3f29713cd
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
# ╠═fc7dbd2e-f621-4953-bd81-b565d74f0af8
# ╠═2f54285c-50cf-49e4-b2b0-cf9aa5fdf585
# ╟─25c85fd1-1921-4707-8185-735751c0d914
# ╟─10236da6-fcef-49bf-9c3c-21f75e83f554
# ╟─adacf377-bbde-4367-b48c-6c4a4721dbfe
# ╠═8f635abb-c9ea-45ca-81ae-e679e0ba923d
# ╠═9b90ddc4-cee5-4063-9f83-ecbff8c7596c
# ╠═7b8bbd0f-c650-49ef-b8d9-c42cd6d8da9e
# ╠═8363c723-74e0-49bf-88ae-e164919f4681
# ╟─3424b05c-e792-4cb1-91ec-d27a04c8d4a2
# ╠═b93903cb-24c9-4ff1-8e31-ee6cb9039c24
# ╠═20c98896-8ae1-4315-a517-3ecf69ca0123
# ╠═d3a39f1d-879e-4e30-9098-6cd8cc58a5af
# ╠═7aad257a-e6b6-4249-b6b2-2834e5be86d5
# ╠═7771dc15-9365-4722-bffc-7ba2534afaf1
# ╟─dba10e9c-1ff0-41fa-af37-86f02cd47ce2
# ╠═e33780b3-1e4c-40ac-a761-8877c79ac170
# ╠═720da1fc-136d-4693-ba1c-b9e7c885e33e
# ╠═c8213569-f42a-49ba-bcca-34d7d5b2b04d
# ╟─7dd5d141-f6a5-467a-ac59-d766254a0e0d
# ╟─e27b5cce-ce85-4e21-85fe-a8757283282f
# ╟─4ac9913e-5d2f-41b3-b782-4d31b441602e
# ╟─d1469e76-1e73-42f5-afd3-90f1bbfb4dfa
# ╟─53ea2aff-3bcc-41ed-be7e-48b8e974f2ba
