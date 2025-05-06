### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 34fd1713-4d0a-4bc9-81e1-bacf418747a2
using Pkg; Pkg.activate("../..")

# ╔═╡ c379b857-4a95-42c0-92bd-bf9df430e1e8
using PlantModules

# ╔═╡ 32b33d71-281e-4d07-b4c1-86f7be6997bf
using PlantGraphs, MultiScaleTreeGraph

# ╔═╡ 8f5728b1-fa78-492e-beab-9e40ebc633ef
using ModelingToolkit, OrdinaryDiffEq, Unitful

# ╔═╡ dc47be70-0bc4-452d-a128-73ff5fce56ce
using Plots

# ╔═╡ 41b55366-5b6e-4489-b37d-342d57fd7b41
using Memoization

# ╔═╡ 876006eb-b8b6-46dc-8077-8f64991c5616
using PlantBiophysics, PlantBiophysics.PlantMeteo, PlantBiophysics.PlantSimEngine

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
plantXEG = PlantModules.readXEG("./structures/beech.xeg");

# ╔═╡ 16f8f64f-130c-4402-87b4-1dc6e7219928
md"""
For manipulating the graph, working in the MultiScaleTreeGraph.jl format is the easiest.
"""

# ╔═╡ d98e14b1-d636-48d2-8825-793166f093f2
md"""
First off, we need to remove all the turtle commands defined in this graph. These are nodes in the graph used to define how the plant should be visualized, but are not actual structural plant parts. Modeling the plant will be easier after removing them. We can easily discern these nodes by the fact they have no dimensions defined (no diameter, length or width).
"""

# ╔═╡ b53c20bb-f25a-4def-a4ed-1ddd35aa30bf
md"""
Secondly, we need to combine the separate dimensions defined in the XEG file into one vector `D`, as is expected by the pre-defined functional modules we'll use.
"""

# ╔═╡ 3499476b-4764-43d8-9886-5e493fb024c5
md"""
Lastly, let's also rename the "ShortShoot" nodetype to simply "Shoot":
"""

# ╔═╡ 5a3a91d3-c90f-46b4-a614-d253e46fb81b
function combine_dimensions(l, d, w)
	if all(isnothing.([l, d, w]))
		return nothing
	elseif isnothing(w) # no width defined => shape is a cylinder
		return 1e2*[d/2, l] # 1e2 to go from m to cm
	else # otherwise we're dealing with a cuboid
		return 1e2*[l, w, 5e-4] # leaf thickness assumed to be 0.5 mm
	end
end;

# ╔═╡ d18e6205-85a4-4488-9e7a-7dfab43954b2
begin
	mtg = PlantModules.convert_to_MTG(plantXEG);
	plant_graph = delete_nodes!(mtg, filter_fun = node -> all(isnothing.([node.diameter, node.length, node.width])))
	transform!(plant_graph, [:length, :diameter, :width] => combine_dimensions => :D)
	traverse!(plant_graph, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")
end;

# ╔═╡ 22bbed9a-f671-4fa1-89b5-c9a80da2e557
DataFrame(plant_graph, [:D])

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
struct_connections = PlantStructure(graphs, intergraph_connections);

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
get_PAR_flux(t) = max(0, 40 * sin(t/24*2*pi - 8));

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
	m_out = run!(m, meteo)
	return only(m_out[:A]) |> x -> max(x, 0) # extract result of the first (and only) timestep
end;

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

# ╔═╡ ea772832-1f9d-4582-a822-83aa0893fdcd
import PlantModules: t, d

# ╔═╡ acda9c67-9918-4c3d-8cb5-df3cceb6c29b
function photosynthesis_module(; name, T, M, shape)	
	@constants (
		uc1 = (10^-6 * 10^-4 * 60^2), 
		[description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)", unit = u"(mol/cm^2/hr) / (µmol/m^2/s)"],
	)
	@parameters (
		T = T, [description = "Temperature", unit = u"K"],
		LAI = 1.0, [description = "Leaf Area Index", unit = u"cm^2 / cm^2"],
		k = 0.5, [description = "Light extinction coefficient", unit = u"N/N"],
		osmotic_frac = 0.3, [description = "Fraction of assimilated carbon which becomes osmotically active solute", unit = u"mol/mol"],
		carbon_decay_rate = 0.2, [description = "Rate at which carbon is consumed for growth", unit = u"hr^-1"],
	)

	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", 
				   unit = u"mol / cm^3"],
		PF(t), [description = "Incoming PAR flux", unit = u"J / s / m^2"],
		A(t), [description = "Carbon assimilation rate", unit = u"µmol / m^2 / s"],
		D(t)[1:length(shape.ϵ_D)], [description = "Dimensions of compartment", 
									unit = u"cm"],
    )

	leafarea(::Cuboid, D::AbstractArray) = D[1] * D[2]

    eqs = [
		PF ~ get_PAR_flux(t)
		A ~ get_assimilation_rate(PF, T, LAI, k)
        d(M) ~ uc1 * osmotic_frac * A * leafarea(shape, D) / volume(shape, D) - carbon_decay_rate*M
    ]
    return ODESystem(eqs, t; name, checks = false)
end;

# ╔═╡ 4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
md"""
#### Defining parameter and initial values 
"""

# ╔═╡ 6456a7ec-47d4-4032-9dec-f9733caa5f66
md"""
This part is again largely the same as in the previous tutorial, except that we have a new functional module to add default values for!
"""

# ╔═╡ 6ffd0511-639d-4304-af8f-93e277e2a4a8
module_defaults = Dict(
	:Internode => Dict(:shape => Cylinder([0.1, 0.1], [0.002, 0.0001]), :M => 300e-6, :K_s => 200),
	:Shoot => Dict(:shape => Cylinder([0.1, 0.1], [0.002, 0.0001]), :M => 350e-6, :K_s => 50),
	:Leaf => Dict(:shape => Cuboid([1.0, 1.0, 1.0], [0.002, 0.002, 5e-4]), :M => 200e-6, :K_s => 5e-5),
	:Soil => Dict(:W_max => 1e4, :T => 293.15), #! W_max
	:Air => Dict(:K => 1e-1)
);

# ╔═╡ a3551081-5da1-4ed4-b3df-9b8929a46f33
connecting_modules = [
	(:Soil, :Internode) => (const_hydraulic_connection, Dict(:K => 1_000)),
    (:Internode, :Internode) => (hydraulic_connection, Dict()),
	(:Internode, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Leaf) => (const_hydraulic_connection, Dict()),
	(:Internode, :Leaf) => (const_hydraulic_connection, Dict()),
    (:Leaf, :Air) => (hydraulic_connection, Dict())
];

# ╔═╡ 7958cad6-3454-4c07-a6b7-14457a49b162
func_connections = PlantFunctionality(; module_defaults, connecting_modules);

# ╔═╡ 930e7ed8-0bfe-4e5a-8890-a1d1ce155881
md"""
### Coupling functional and structural modules
"""

# ╔═╡ dc8cffc1-675f-43ce-aadd-c8b786bf5889
md"""
The same goes for this section: all we have to change compared to previous tutorial is coupling our new photosynthesis module to Leaf nodes instead of the constant carbon module.
"""

# ╔═╡ 77686779-bb3c-4d9f-b95c-2613a74e4444
module_coupling = Dict(
	:Internode => [hydraulic_module, constant_carbon_module, K_module],
	:Shoot => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, photosynthesis_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
);

# ╔═╡ 210d81ef-153e-4744-8266-80af4099770c
md"""
### Bringing it all together
"""

# ╔═╡ bc7573e7-bcd6-4347-9b0c-9111a824c9b5
md"""

"""

# ╔═╡ 7090562e-92e1-409a-9a91-856366e3a20e
system = generate_system(struct_connections, func_connections, module_coupling,
		checkunits = false);

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
time_span = (0, 7*24.0); # We'll simulate our problem for a timespan of 7 days

# ╔═╡ 50d6fc31-80f5-4db7-b716-b26765008a0d
prob = ODEProblem(sys_simpl, [], (0.0, 5*24));

# ╔═╡ c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
sol = solve(prob);

# ╔═╡ a6608eff-9399-443c-a33a-c62341f7b14c
md"""
## Answering the toy problem

Just like last tutorial, we could just plot the soil's relative water content. However, there's a lot of interesting other variables we could look at, like those relating to the carbon dynamics we've just introduced!
"""

# ╔═╡ c616121f-ece2-4a3e-9743-32ef70b1c493
plotgraph(sol, graphs[1], varname = :W)

# ╔═╡ c65f22fd-82b2-462c-b350-8906495963c6
plotgraph(sol, graphs[1:2], varname = :Ψ)

# ╔═╡ 52aa2fee-aaec-4bbf-b753-eeafb93809b3
plotgraph(sol, graphs[2], varname = :W_r)

# ╔═╡ d48f25b1-8e6b-42aa-a2c3-be3d85b3b8d2
md"""
There's a lot more water in the soil now, so the relative water content goes down a lot slower.
"""

# ╔═╡ fd5e242f-25b6-46e8-bb6f-4112073c2bc4
plotgraph(sol, plant_graph, varname = :M, structmod = :Leaf)

# ╔═╡ 2009c5fc-f237-4bf3-af5a-c0adee2c5aed
md"""
We defined the same input radiation, assimilation rate and thickness for all leaves, so it's no surprise they all have a similar metabolite concentration.
"""

# ╔═╡ fe4df2d4-878e-41aa-8860-991c891e2dd2
plotgraph(sol, plant_graph, varname = :W, structmod = :Leaf)

# ╔═╡ 2470d664-ccc9-46ac-9964-b39dabe2ce1b
md"""
We can see that the water content over time follows a periodic pattern, which is what we would expect.
"""

# ╔═╡ 9aea25dc-8543-4249-aada-df02fe93527e
plotgraph(sol, plant_graph, varname = :W, structmod = :Internode)

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
# ╠═dc47be70-0bc4-452d-a128-73ff5fce56ce
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╟─aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
# ╟─6ef5c63a-b753-43ae-baee-f6c24313a385
# ╠═e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
# ╟─16f8f64f-130c-4402-87b4-1dc6e7219928
# ╟─d98e14b1-d636-48d2-8825-793166f093f2
# ╟─b53c20bb-f25a-4def-a4ed-1ddd35aa30bf
# ╟─3499476b-4764-43d8-9886-5e493fb024c5
# ╠═5a3a91d3-c90f-46b4-a614-d253e46fb81b
# ╠═d18e6205-85a4-4488-9e7a-7dfab43954b2
# ╠═22bbed9a-f671-4fa1-89b5-c9a80da2e557
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
# ╠═41b55366-5b6e-4489-b37d-342d57fd7b41
# ╠═876006eb-b8b6-46dc-8077-8f64991c5616
# ╠═f5a529ff-90d9-4601-97df-66c24a7faaa2
# ╟─899b48f7-aed8-4f20-a139-35196940539c
# ╠═a417f515-6cad-486a-8ce0-c1eb72bb5939
# ╟─2508d5d8-d4a2-440a-94e4-b17d6a5dc4dd
# ╠═ea772832-1f9d-4582-a822-83aa0893fdcd
# ╠═acda9c67-9918-4c3d-8cb5-df3cceb6c29b
# ╟─4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
# ╟─6456a7ec-47d4-4032-9dec-f9733caa5f66
# ╠═6ffd0511-639d-4304-af8f-93e277e2a4a8
# ╠═a3551081-5da1-4ed4-b3df-9b8929a46f33
# ╠═7958cad6-3454-4c07-a6b7-14457a49b162
# ╟─930e7ed8-0bfe-4e5a-8890-a1d1ce155881
# ╟─dc8cffc1-675f-43ce-aadd-c8b786bf5889
# ╠═77686779-bb3c-4d9f-b95c-2613a74e4444
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
# ╠═c616121f-ece2-4a3e-9743-32ef70b1c493
# ╠═c65f22fd-82b2-462c-b350-8906495963c6
# ╠═52aa2fee-aaec-4bbf-b753-eeafb93809b3
# ╟─d48f25b1-8e6b-42aa-a2c3-be3d85b3b8d2
# ╠═fd5e242f-25b6-46e8-bb6f-4112073c2bc4
# ╟─2009c5fc-f237-4bf3-af5a-c0adee2c5aed
# ╠═fe4df2d4-878e-41aa-8860-991c891e2dd2
# ╟─2470d664-ccc9-46ac-9964-b39dabe2ce1b
# ╠═9aea25dc-8543-4249-aada-df02fe93527e
# ╟─91d295df-e76d-402b-8d59-68a7588b57ae
