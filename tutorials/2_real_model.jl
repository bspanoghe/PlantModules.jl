### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 34fd1713-4d0a-4bc9-81e1-bacf418747a2
using Pkg; Pkg.activate(".")

# ╔═╡ c379b857-4a95-42c0-92bd-bf9df430e1e8
using PlantModules

# ╔═╡ 32b33d71-281e-4d07-b4c1-86f7be6997bf
using PlantGraphs, MultiScaleTreeGraph

# ╔═╡ 8f5728b1-fa78-492e-beab-9e40ebc633ef
using ModelingToolkit, DifferentialEquations, Unitful

# ╔═╡ 876006eb-b8b6-46dc-8077-8f64991c5616
using PlantBiophysics, PlantBiophysics.PlantMeteo, PlantBiophysics.PlantSimEngine

# ╔═╡ 41b55366-5b6e-4489-b37d-342d57fd7b41
using Memoization

# ╔═╡ dc47be70-0bc4-452d-a128-73ff5fce56ce
using Plots; import GLMakie.draw

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 2: Making a real model
"""

# ╔═╡ 6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
md"""
## Introduction 👋

Last tutorial has covered the most basic functionality of the package in order to get a model running of the hydraulics in a very simple plant. In this tutorial, we'll expand on the same idea but cover some additional functionality of PlantModules to make a more useful model. Specifically, the following will be discussed:
- Working with larger plant structures
- Adding new functional processes
- Integrating observational data into the model
"""

# ╔═╡ 1144887a-a4c7-46f6-9cf8-cba50cf873d0
md"""
### Toy problem description

For the second tutorial, we'll again keep track of the water dynamics between the soil and a plant. However, this time we're scaling up: instead of a pepper seedling, we'll be considering a beech tree. Additionally, we'll include photosynthesis into the model this time combined with (simulated) weather data!

![plantfigu](https://www.woodlandtrust.org.uk/media/4550/copper-beech-tree-mature-alamy-e08fat-andrew-roland.jpg)
"""
#! replace picture with something royalty free. a drawing?

# ╔═╡ 62afb755-28a3-46c0-be89-4766a46d789a
md"""
(Note: The plant we'll be simulating is a little smaller than the one in the picture above for the sake of reducing runtime)
"""

# ╔═╡ cbbb6235-4482-4aec-85bf-3f315e8b9df7
md"""
## Loading packages
"""

# ╔═╡ 0255cc11-11a2-4c72-b676-26b2f0ace488
import .PlantModules: Cilinder, Cuboid, volume, convert_to_PG, convert_to_MTG, generate_system #! Are exported, why don't they work >:(

# ╔═╡ b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
md"""
## Defining the structural modules
"""

# ╔═╡ aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
md"""
### Working with real plant structures
"""

# ╔═╡ 6ef5c63a-b753-43ae-baee-f6c24313a385
md"""
For practical applications, it is generally not adviced to write out the entire plant structure by hand. An easier way to acquire them is by use of L-systems or some other kind of rewriting system. This can be done in Julia, making use of e.g. [PlantGeom.jl](https://github.com/VEZY/PlantGeom.jl), or using other software, in which case the plant structure needs to be saved as a file. For this tutorial we'll consider the latter case and use a Beech structure generated in [GroIMP](https://wwwuser.gwdguser.de/~groimp/grogra.de/software/groimp/index.html) from its [example library](https://groimp.wordpress.com/fspm/). The plant structure was saved as an XEG file, a file format for FSPMs created by GroIMP and OpenAlea. Another popular file format for plant structures is the MTG (Multiscale Tree Graph) file, which can be read using the MultiScaleTreeGraph.jl package.
"""

# ╔═╡ e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
plantXEG = PlantModules.readXEG("./structures/beech.xeg")

# ╔═╡ 1d3d04e4-7bb2-4dba-a25a-5bcde565ee62
md"""
We can then convert the graph to other formats.

`convert_to_PG` converts a graph to PlantGeom.jl format, which can be used for its visualisation functions, for example:
"""

# ╔═╡ d15c0747-6e2a-4b61-9b83-12cc4585cbf9
PlantModules.convert_to_PG(plantXEG) |> draw

# ╔═╡ 16f8f64f-130c-4402-87b4-1dc6e7219928
md"""
For manipulating the graph, working in the MultiScaleTreeGraph.jl format is the easiest.
"""

# ╔═╡ 69480c82-62c9-48d5-bb9e-3d698603c937
mtg = PlantModules.convert_to_MTG(plantXEG);

# ╔═╡ b1c4c7f6-8a54-41f6-ae08-7521f5b9333a
DataFrame(mtg, [:diameter, :length, :width])

# ╔═╡ d98e14b1-d636-48d2-8825-793166f093f2
md"""
First off, we need to remove all the turtle commands defined in this graph. These are nodes in the graph used to define how the plant should be visualized, but are not actual structural plant parts. Modeling the plant will be easier after removing them. We can easily discern these nodes by the fact they have no dimensions defined (no diameter, length or width).
"""

# ╔═╡ 0065bb33-2fb8-465e-b908-ba0df2248e18
plant_graph = delete_nodes!(mtg, filter_fun = node -> all(isnothing.([node.diameter, node.length, node.width])))

# ╔═╡ 1b87aab5-6fb5-446b-8acd-e64d9c00754b
DataFrame(plant_graph, [:diameter, :length, :width])

# ╔═╡ b53c20bb-f25a-4def-a4ed-1ddd35aa30bf
md"""
Secondly, we need to combine the separate dimensions defined in the XEG file into one vector `D`, as is expected by the pre-defined functional modules we'll use.
"""

# ╔═╡ 5a3a91d3-c90f-46b4-a614-d253e46fb81b
function combine_dimensions(l, d, w)
	if isnothing(w) # no width defined => plant part is a cilinder
		return [100*d/2, 100*l] # file has dimensions in m, our functional modules expect cm
	else # otherwise plant part is a cuboid
		return [100*l, 100*w, 0.03] # leaf thickness not defined => set to 0.03 cm
	end
end

# ╔═╡ ffd91a3f-1475-471d-be28-98f5aaaacabb
transform!(plant_graph, [:length, :diameter, :width] => combine_dimensions => :D)

# ╔═╡ 3499476b-4764-43d8-9886-5e493fb024c5
md"""
Lastly, let's also rename the "ShortShoot" nodetype to simply "Shoot":
"""

# ╔═╡ 2919dad7-f954-40b3-a47b-20ae62cdce1b
traverse!(plant_graph, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")

# ╔═╡ 22bbed9a-f671-4fa1-89b5-c9a80da2e557
DataFrame(plant_graph, [:D])

# ╔═╡ 068741d7-7d28-4376-9575-2a0f29675928
convert_to_PG(plant_graph) |> draw

# ╔═╡ 98eac4c4-b39a-4e11-917a-90b03d7385d1
md"""
#### The environment

For this part, we'll just do the same as last tutorial.
"""

# ╔═╡ e00c5135-1d66-4dec-8283-40ebe06a8038
struct Soil <: PlantGraphs.Node end

# ╔═╡ dac02191-b640-40f5-a7d6-e6b06b946c23
struct Air <: PlantGraphs.Node end

# ╔═╡ e7976ec8-c654-4ee9-a01d-1f961a12c0c3
intergraph_connections = [[1, 2] => (plant_graph, :Soil), [1, 3] => (:Leaf, :Air)];

# ╔═╡ ef827a9c-6447-410e-8b0d-320ba16c4137
graphs = [plant_graph, Soil(), Air()];

# ╔═╡ d5c92ba1-0d5a-4d3c-8dd6-c760e6d2a67c
struct_connections = [graphs, intergraph_connections];

# ╔═╡ 43211f69-6bfe-4fd1-b474-65d0601558de
md"""
### Defining the functional modules

Previous tutorial, we simply made use of the pre-implemented functional modules. Now, let's take a look at how we can define our own. We'll start by making a simple model for carbon dynamics.
"""

# ╔═╡ e4b9805e-0d4e-48da-a80c-d5cff53c3f15
md"""
First thing we need is to define how much photosynthetically active radiation (PAR) our leaves receive. We'd like to define this as a general Julia function rather than a differential equation so that we have more freedom in our definition. This allows you to, for example, write a function which interpolates between observed PAR values. We'll just use a sine function bound to the positive values to simulate a simple day-night cycle:
"""

# ╔═╡ 56a3ccbc-8d31-4579-ad36-15882072ebd6
get_PAR_flux(t) = max(0, 400 * sin(t/24*2*pi - 8))

# ╔═╡ 00a71441-548d-4172-91a8-6169285a1b11
@register_symbolic get_PAR_flux(t)

# ╔═╡ aeac60ec-1760-4f21-b68e-02c747356af3
md"""
The second ingredient we'd like to define outside of our differential equations is the actual photosynthesis model. We _could_ define this ourselves using differential equations, but we could also simply use an implementation already available in the Julia ecosystem, namely from the wonderful [PlantBiophysics.jl](https://github.com/VEZY/PlantBiophysics.jl). All we have to do is wrap this model in a julia function with the desired inputs and outputs and we're done.
"""

# ╔═╡ f5a529ff-90d9-4601-97df-66c24a7faaa2
@memoize function get_assimilation_rate(PAR_flux, T, LAI, k)
	Kelvin_to_C = -273.15
	meteo = Atmosphere(T = T + Kelvin_to_C, Wind = 1.0, P = 101.3, Rh = 0.65, Ri_PAR_f = PAR_flux)
	m = ModelList(
		Fvcb(), # calculate CO2 assimilation rate
		Medlyn(0.03, 0.92), # calculate stomatal conductance, see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2486.2010.02375.x
		Beer(k), # calculate amount of light intercepted
		status = (Tₗ = meteo[:T], LAI = LAI, Cₛ = meteo[:Cₐ], Dₗ = meteo[:VPD], RI_PAR_f = meteo[:Ri_PAR_f])
	)
	run!(m, meteo)
	return only(m[:A]) |> x -> max(x, 0) # extract result of the first (and only) timestep
end

# ╔═╡ 899b48f7-aed8-4f20-a139-35196940539c
md"""
Note the `@memoize` in front of the function. This macro comes from `Memoization.jl`, and is included here because it speeds up the model by a large factor. Simply put, `@memoize` will make the function remember which inputs it has seen before and what the corresponding outputs are. When it encounters a previously seen set of inputs, it will then simply return these outputs without running the entire function again. This can save a lot of time for costly functions such as this photosynthesis model.
"""

# ╔═╡ a417f515-6cad-486a-8ce0-c1eb72bb5939
@register_symbolic get_assimilation_rate(PAR_flux, T, LAI, k)

# ╔═╡ 2508d5d8-d4a2-440a-94e4-b17d6a5dc4dd
md"""
Now for to define the ModelingToolkit function!
"""

# ╔═╡ acda9c67-9918-4c3d-8cb5-df3cceb6c29b
function photosynthesis_module(; name, T, M, shape)
	@variables t, [description = "Time", unit = u"hr"]
	d = Differential(t)
	
	@constants (
		uc1 = (10^-6 * 10^-4 * 60^2), 
		[description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)", unit = u"(mol/cm^2/hr) / (µmol/m^2/s)"],
	)
	@parameters (
		T = T, [description = "Temperature", unit = u"K"],
		LAI = 2.0, [description = "Leaf Area Index", unit = u"cm^2 / cm^2"],
		k = 0.5, [description = "Light extinction coefficient", unit = u"N/N"],
		osmotic_frac = 0.3, [description = "Fraction of assimilated carbon which becomes osmotically active solute", unit = u"mol/mol"],
		carbon_decay_rate = 0.1, [description = "Rate at which carbon is consumed for growth", unit = u"hr^-1"],
	)

	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / cm^3"],
		PF(t), [description = "Incoming PAR flux", unit = u"J / s / m^2"],
		A(t), [description = "Carbon assimilation rate", unit = u"µmol / m^2 / s"],
		D(t)[1:length(shape.ϵ_D)], [description = "Dimensions of compartment", unit = u"cm"],
    )

	leafarea(::Cuboid, D::AbstractArray) = D[1] * D[2]

    eqs = [
		PF ~ get_PAR_flux(t)
		A ~ get_assimilation_rate(PF, T, LAI, k)
        d(M) ~ uc1 * osmotic_frac * A * leafarea(shape, D) / volume(shape, D) - carbon_decay_rate*M
    ]
    return ODESystem(eqs, t; name, checks = false)
end

# ╔═╡ 4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
md"""
#### Defining parameter and initial values 
"""

# ╔═╡ 6456a7ec-47d4-4032-9dec-f9733caa5f66
md"""
This part is again largely the same as in the previous tutorial, except that we have a new functional module to add default values for!
"""

# ╔═╡ d67a31c1-5325-4112-8e3f-832a317f1594
begin

C_stem = 300e-6
C_shoot = 350e-6
C_leaf = 400e-6

default_params = merge(PlantModules.default_params, 
	(photosynthesis_module = (shape = Cuboid(ϵ_D = [0, 0, 0], ϕ_D = [0, 0, 0]), T = 293.15),),
	(waterdependent_hydraulic_connection = (K_max = 0,),),
)

default_u0s = merge(PlantModules.default_u0s,
	(photosynthesis_module = (M = 0,),),
	(waterdependent_hydraulic_connection = (),)
)

module_defaults = (
	Internode = (shape = Cilinder(ϵ_D = [10.0, 300.0], ϕ_D = [1e-3, 1e-5]), M = C_stem),
	Shoot = (shape = Cilinder(ϵ_D = [10.0, 300.0], ϕ_D = [3e-5, 3e-2]), M = C_shoot),
	Leaf = (shape = Cuboid(ϵ_D = [15.0, 10.0, 1000.0], ϕ_D = [3e-3, 3e-3, 1e-5]), M = C_leaf),
	Soil = (W_max = 5000.0, T = 293.15, W_r = 0.9),
	Air = (W_r = 0.8,)
)
	
end

# ╔═╡ 930e7ed8-0bfe-4e5a-8890-a1d1ce155881
md"""
### Coupling functional and structural modules
"""

# ╔═╡ dc8cffc1-675f-43ce-aadd-c8b786bf5889
md"""
The same goes for this section: all we have to change compared to previous tutorial is coupling our new photosynthesis module to Leaf nodes instead of the constant carbon module.
"""

# ╔═╡ 77686779-bb3c-4d9f-b95c-2613a74e4444
module_coupling = [
	PlantModules.hydraulic_module => [:Internode, :Shoot, :Leaf],
	PlantModules.constant_carbon_module => [:Internode, :Shoot],
	photosynthesis_module => [:Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air],
]

# ╔═╡ 4755a0a2-f585-4000-8a55-f4512a582281
md"""
### Defining the connections
"""

# ╔═╡ 95f47c11-2dfe-41bd-b24b-4e1a0bd780e9
md"""
Since we immediately read in the structure of the plant from a file, the structural connections are already defined. All that remains is the define the functional information of the edges.
"""

# ╔═╡ 87821d4e-253c-4d63-ab73-154001e9879d
connecting_modules = [
	(:Soil, :Internode) => (PlantModules.hydraulic_connection, [:K => 15]),
    (:Internode, :Internode) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Internode, :Shoot) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Shoot, :Shoot) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Shoot, :Leaf) => (PlantModules.hydraulic_connection, [:K => 3]),
	(:Internode, :Leaf) => (PlantModules.hydraulic_connection, [:K => 3]),
    (:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 3e-3])
]

# ╔═╡ 9be1023a-d117-43e0-9ba6-bc9f0589024e
func_connections = [connecting_modules, PlantModules.multi_connection_eqs]

# ╔═╡ 210d81ef-153e-4744-8266-80af4099770c
md"""
### Bringing it all together
"""

# ╔═╡ bc7573e7-bcd6-4347-9b0c-9111a824c9b5
md"""

"""

# ╔═╡ 7090562e-92e1-409a-9a91-856366e3a20e
system = generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
);

# ╔═╡ f14d7782-326b-4b1e-a35a-ea479f76e85c
sys_simpl = structural_simplify(system);

# ╔═╡ d51795b2-32d3-455c-b785-5e501cfbdd08
md"""
Calling this constructor will create a `PlantSystem`-type variable containing all required information for running the model with ModelingToolkit. It is possible to fine-tune the model even further at this stage as described in the [Customizing the model](nothinghere) section of the docs, thought this should generally not be required.
"""

# ╔═╡ d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
md"""
## Running the model 🏃‍♂️
"""

# ╔═╡ 6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
md"""
The rest of the modeling workflow is mostly taken care of by the ModelingToolkit and DifferentialEquations Julia packages, with some syntactic sugar added by PlantModules. For users that are unfamiliar with the package, it is recommended to take a brief look at [the ModelingToolkit docs](https://docs.sciml.ai/ModelingToolkit/stable/) before proceeding.
"""

# ╔═╡ bf114636-1e35-49f1-9407-f472b443a9ea
time_span = (0, 7*24.0) # We'll simulate our problem for a timespan of 7 days

# ╔═╡ 50d6fc31-80f5-4db7-b716-b26765008a0d
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), time_span)

# ╔═╡ c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
@time sol = solve(prob);

# ╔═╡ a6608eff-9399-443c-a33a-c62341f7b14c
md"""
## Answering the toy problem

Just like last tutorial, we could just plot the soil's relative water content. However, there's a lot of interesting other variables we could look at, like those relating to the carbon dynamics we've just introduced!
"""

# ╔═╡ 52aa2fee-aaec-4bbf-b753-eeafb93809b3
PlantModules.plotgraph(sol, graphs[2], func_varname = :W_r)

# ╔═╡ d48f25b1-8e6b-42aa-a2c3-be3d85b3b8d2
md"""
There's a lot more water in the soil now, so the relative water content goes down a lot slower.
"""

# ╔═╡ fd5e242f-25b6-46e8-bb6f-4112073c2bc4
PlantModules.plotgraph(sol, plant_graph, func_varname = :M, struct_module = :Leaf)

# ╔═╡ 2009c5fc-f237-4bf3-af5a-c0adee2c5aed
md"""
We defined the same input radiation, assimilation rate and thickness for all leaves, so it's no surprise they all have the same metabolite concentration.
"""

# ╔═╡ fe4df2d4-878e-41aa-8860-991c891e2dd2
PlantModules.plotgraph(sol, plant_graph, func_varname = :W, struct_module = :Leaf)

# ╔═╡ 2470d664-ccc9-46ac-9964-b39dabe2ce1b
md"""
We can see that the water content over time follows a periodic pattern, which is what we would expect.
"""

# ╔═╡ 9aea25dc-8543-4249-aada-df02fe93527e
PlantModules.plotgraph(sol, plant_graph, func_varname = :W, struct_module = :Internode)

# ╔═╡ 91d295df-e76d-402b-8d59-68a7588b57ae
md"""
The same goes for the internodes, since their water exchange rate with the leaves are also based on those metabolite concentrations through the osmotic potential.
"""

# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
# ╟─1144887a-a4c7-46f6-9cf8-cba50cf873d0
# ╟─62afb755-28a3-46c0-be89-4766a46d789a
# ╟─cbbb6235-4482-4aec-85bf-3f315e8b9df7
# ╠═34fd1713-4d0a-4bc9-81e1-bacf418747a2
# ╠═c379b857-4a95-42c0-92bd-bf9df430e1e8
# ╠═32b33d71-281e-4d07-b4c1-86f7be6997bf
# ╠═8f5728b1-fa78-492e-beab-9e40ebc633ef
# ╠═876006eb-b8b6-46dc-8077-8f64991c5616
# ╠═41b55366-5b6e-4489-b37d-342d57fd7b41
# ╠═dc47be70-0bc4-452d-a128-73ff5fce56ce
# ╠═0255cc11-11a2-4c72-b676-26b2f0ace488
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╟─aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
# ╟─6ef5c63a-b753-43ae-baee-f6c24313a385
# ╠═e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
# ╟─1d3d04e4-7bb2-4dba-a25a-5bcde565ee62
# ╠═d15c0747-6e2a-4b61-9b83-12cc4585cbf9
# ╟─16f8f64f-130c-4402-87b4-1dc6e7219928
# ╠═69480c82-62c9-48d5-bb9e-3d698603c937
# ╠═b1c4c7f6-8a54-41f6-ae08-7521f5b9333a
# ╟─d98e14b1-d636-48d2-8825-793166f093f2
# ╠═0065bb33-2fb8-465e-b908-ba0df2248e18
# ╠═1b87aab5-6fb5-446b-8acd-e64d9c00754b
# ╟─b53c20bb-f25a-4def-a4ed-1ddd35aa30bf
# ╠═5a3a91d3-c90f-46b4-a614-d253e46fb81b
# ╠═ffd91a3f-1475-471d-be28-98f5aaaacabb
# ╟─3499476b-4764-43d8-9886-5e493fb024c5
# ╠═2919dad7-f954-40b3-a47b-20ae62cdce1b
# ╠═22bbed9a-f671-4fa1-89b5-c9a80da2e557
# ╠═068741d7-7d28-4376-9575-2a0f29675928
# ╟─98eac4c4-b39a-4e11-917a-90b03d7385d1
# ╠═e00c5135-1d66-4dec-8283-40ebe06a8038
# ╠═dac02191-b640-40f5-a7d6-e6b06b946c23
# ╠═e7976ec8-c654-4ee9-a01d-1f961a12c0c3
# ╠═ef827a9c-6447-410e-8b0d-320ba16c4137
# ╠═d5c92ba1-0d5a-4d3c-8dd6-c760e6d2a67c
# ╟─43211f69-6bfe-4fd1-b474-65d0601558de
# ╟─e4b9805e-0d4e-48da-a80c-d5cff53c3f15
# ╠═56a3ccbc-8d31-4579-ad36-15882072ebd6
# ╠═00a71441-548d-4172-91a8-6169285a1b11
# ╟─aeac60ec-1760-4f21-b68e-02c747356af3
# ╠═f5a529ff-90d9-4601-97df-66c24a7faaa2
# ╟─899b48f7-aed8-4f20-a139-35196940539c
# ╠═a417f515-6cad-486a-8ce0-c1eb72bb5939
# ╟─2508d5d8-d4a2-440a-94e4-b17d6a5dc4dd
# ╠═acda9c67-9918-4c3d-8cb5-df3cceb6c29b
# ╟─4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
# ╟─6456a7ec-47d4-4032-9dec-f9733caa5f66
# ╠═d67a31c1-5325-4112-8e3f-832a317f1594
# ╟─930e7ed8-0bfe-4e5a-8890-a1d1ce155881
# ╟─dc8cffc1-675f-43ce-aadd-c8b786bf5889
# ╠═77686779-bb3c-4d9f-b95c-2613a74e4444
# ╟─4755a0a2-f585-4000-8a55-f4512a582281
# ╟─95f47c11-2dfe-41bd-b24b-4e1a0bd780e9
# ╠═87821d4e-253c-4d63-ab73-154001e9879d
# ╠═9be1023a-d117-43e0-9ba6-bc9f0589024e
# ╟─210d81ef-153e-4744-8266-80af4099770c
# ╟─bc7573e7-bcd6-4347-9b0c-9111a824c9b5
# ╠═7090562e-92e1-409a-9a91-856366e3a20e
# ╠═f14d7782-326b-4b1e-a35a-ea479f76e85c
# ╟─d51795b2-32d3-455c-b785-5e501cfbdd08
# ╟─d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
# ╟─6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
# ╠═bf114636-1e35-49f1-9407-f472b443a9ea
# ╠═50d6fc31-80f5-4db7-b716-b26765008a0d
# ╠═c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
# ╟─a6608eff-9399-443c-a33a-c62341f7b14c
# ╠═52aa2fee-aaec-4bbf-b753-eeafb93809b3
# ╟─d48f25b1-8e6b-42aa-a2c3-be3d85b3b8d2
# ╠═fd5e242f-25b6-46e8-bb6f-4112073c2bc4
# ╟─2009c5fc-f237-4bf3-af5a-c0adee2c5aed
# ╠═fe4df2d4-878e-41aa-8860-991c891e2dd2
# ╟─2470d664-ccc9-46ac-9964-b39dabe2ce1b
# ╠═9aea25dc-8543-4249-aada-df02fe93527e
# ╟─91d295df-e76d-402b-8d59-68a7588b57ae
