### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 57b8dcb8-9baa-4ddf-9368-431b1be5850e
# ╠═╡ show_logs = false
using Pkg; Pkg.activate("../..")

# ╔═╡ 662c0476-70aa-4a60-a81c-d7db2248b728
# ╠═╡ show_logs = false
using PlantModules

# ╔═╡ 65f88593-1180-447a-900f-49aef4647cd1
using PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Plots

# ╔═╡ 7ff01837-8d59-4df7-ad50-08157224e166
using PlutoUI; TableOfContents()

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 1: Package basics
"""

# ╔═╡ 6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
md"""
In this tutorial, we'll cover the basics of this package. The main topics are:
- Introduction: the concept of the package.
- Modeling context
- Creating a modular plant model
  - Defining structural and functional plant modules
  - Specifying model parameters
  - Coupling structure and function
- Running the model
- Visualizing the results
"""

# ╔═╡ 06004b15-17d4-4f43-a6f4-b67f5610ebac
md"## Setup"

# ╔═╡ bc44fae7-6c4a-4f27-be51-dbabf4037601
md"""
Before we can get started, we need to do some basic setup: activating our project and loading in the required packages.
"""

# ╔═╡ aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
md"""
## Introduction: modeling with PlantModules.jl
"""

# ╔═╡ 6ef5c63a-b753-43ae-baee-f6c24313a385
md"""
As the name suggests, one of the most important goals of PlantModules is to enable the user to easily model plant growth in a _modular_ manner. In order to achieve this, we will define our plant in function of a few sets of similarly behaving parts or "modules" for short. Parts of the same module can be similar on the structural - and on the functional level.

On the structural level, some examples are:
- An oak tree can be considered a repeating module in a forest.
- A branch can be considered a repeating module in a tree.
- A collenchyma cell can be considered a repeating module in a branch.

On the functional level, on the other hand:
- Most structural modules of a tree (branches, leaves, root segments, etc.) contain water which will flow between them as dictated by hydraulic laws, which plays an important role in plant growth. All these structural modules share this same functional module. 
- One or more of the plant's structural modules will assimilate carbon through photosynthesis, while the others do not. On the cellular level, only cells containing chloroplasts should share this functional module.

As such, the workflow to create a model in PlantModules boils down to defining these modules and how they interact.
"""

# ╔═╡ 1144887a-a4c7-46f6-9cf8-cba50cf873d0
md"""
## Example problem definition

For the first tutorial, we will model the response of a hypothetical tomato plant to different transpiration rates. In higher plants, an increase in transpiration rate causes a reduction in the water potential of the stem. This phenomenon can be explained by the cohesion-tension theory, which states that the water in plant vessels forms a continuous network through its strong cohesive forces, and pressure reductions caused by leaf transpiration cause pressure gradients throughout the entire network. We will verify here that our basic functional modules can capture this phenomenon.
"""

# ╔═╡ d131a541-9357-4d1a-b6a6-539bf5c9d884
md"## Defining the model"

# ╔═╡ b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
md"""
### Structure
"""

# ╔═╡ aec7bcd6-6f27-4cf5-a955-f4d59e778fd3
md"""
#### The plant

The first step in modeling our plant's growth is defining its structure.
For the first tutorial, we'll consider a very simple example: a tine tomato plant with a small stem and two leaves.
"""

# ╔═╡ 659e911c-8af2-4a66-855a-e333c41120c1
md"""
PlantModules uses **graphs** to define the relationships between modules. As such, structural modules need to be implemented as graph nodes of some sort. We'll be using the graph implementation from [the PlantGraphs](https://virtualplantlab.com/stable/manual/Graphs/) package for this tutorial, though [any graph implementation can be used](nothinghere).
"""

# ╔═╡ 0cc02e82-4fe8-4f27-a2d2-eb4bfba6b291
md"""
For our example, we will define two structural modules: a `Stem` module and a `Leaf` module. We will give both modules an attribute `D`, which corresponds to the name of the variable expressing the dimensions of the structural module. This is a vector containing one to three numbers, depending on whether the plant part is modelled as a sphere, a cylinder or a cuboid. The Leaves will be represented as flat cuboids, while the stem will be represented as a cylinder.
"""

# ╔═╡ 6b7ebc68-f4a1-4ed6-b12b-e4ac5ee9b00a
struct Stem <: Node end

# ╔═╡ d57718c2-c77d-42a8-924f-ebdfcd51d919
struct Leaf <: Node
	D::Vector
end

# ╔═╡ c4dc4961-ba2b-4b23-b80e-7d4eb8d9a9f4
md"We now define the plant structure using a graph. By specifying the value of the attribute `D`, we can specify different initial sizes for the different plant parts. Note that the cylindrical stem only needs 2 values specified (radius and length), while the cuboid leaves need 3 values (length, width and thickness)."

# ╔═╡ 9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
plant_graph = Stem() + Stem() +
	(Leaf([5.0, 3.0, 0.05]), Leaf([1.0, 0.5, 0.05]));

# ╔═╡ 27160915-5943-4048-b363-d123846b02b8
plotstructure(plant_graph)

# ╔═╡ 98eac4c4-b39a-4e11-917a-90b03d7385d1
md"""
#### The environment

Plants need an environment to grow in. For most plants, the most basic environmental compartments that need defining are the soil, from which the plant gets water and nutrients, and the air, with which the plant exchanges gasses. We can make these compartments as complex as we want, such as different soil compartments for different layers, but we will stick to one graph node for each for now.
"""

# ╔═╡ e00c5135-1d66-4dec-8283-40ebe06a8038
struct Soil <: Node end

# ╔═╡ dac02191-b640-40f5-a7d6-e6b06b946c23
struct Air <: Node end

# ╔═╡ 3bf42137-1551-44d6-b7af-eab13a97b6ef
soil_graph = Soil();

# ╔═╡ db8fe96d-c8c2-47c0-9377-281ce577a920
air_graph = Air();

# ╔═╡ c36668fa-3046-4967-a616-841a048ea7f9
md"""
We also need to define the connections between these environmental modules and the plant. Sadly, connecting both leaf nodes to the air would introduce a cycle into our graph, which is something many plant-focused graph implementations do not support.

To overcome this, PlantModules allows the user to input multiple graphs and then specify how to connect them:
"""

# ╔═╡ 8980f1eb-e461-4624-9cce-83a7cf822349
graphs = [plant_graph, soil_graph, air_graph];

# ╔═╡ 61bf737a-2226-42dc-b93a-a8f578048268
intergraph_connections = [[1, 2] => (getnodes(plant_graph)[1], :Soil), [1, 3] => (:Leaf, :Air)];

# ╔═╡ 20049311-d6e6-41d3-a0d8-8bad88c174f9
md"""
The above syntax means: "connect the first node from `graphs[1]` to all Soil nodes from `graphs[2]`", and all `Leaf` nodes from `graphs[1]` to all `Air` nodes from `graphs[3]`. Note that the order matters! For example, `[1, 2] => (:Root, :Soil)` is not the same as `[2, 1] => (:Root, :Soil)`.
"""

# ╔═╡ 668e02ee-ac78-4b3d-983d-402ec04584ef
md"""
The full structural connection information can then be acquired by bundling the graphs and connections into the `PlantStructure` wrapper.
"""

# ╔═╡ caab574a-05c5-4c0d-9ae4-19fd514a6c6c
plantstructure = PlantStructure(graphs, intergraph_connections);

# ╔═╡ 579c41a8-e054-474a-985e-314c9a7657ce
plotstructure(plantstructure)

# ╔═╡ f03a61ce-a0ff-43ff-abdd-2342f76e8c93
md"""
!!! note
	For our tiny example plant, we only have one structural module that actually repeats, somewhat defeating the purpose of modelling in a modular manner: we may as well write this entire model out by hand! Rest assured, however, that the approach we're seeing here also works for larger plants, as we will see in the next tutorial.
"""

# ╔═╡ 43211f69-6bfe-4fd1-b474-65d0601558de
md"""
### Function
"""

# ╔═╡ a6552c16-b013-43af-9483-639a79944749
md"""
Now that we have an idea of the plant's structural modules, we need to assign them some kind of functionality. `PlantModules` expresses these as sets of differential equations, implemented in [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/). Assigning functionality then comes down to writing out the desired functional processes as `ModelingToolkit` models (or choosing them out from the prebuilt ones), mapping them to the structural modules, and finally specifying the initial- and parameter values for the different components of the system.
"""

# ╔═╡ 52bb8235-abea-45c5-929c-5824bc8681c4
md"#### Defining functional modules"

# ╔═╡ c04564c4-4fb5-47bf-bc14-77aaebdece15
md"""
There are two types of functional module: next to the typical functional modules describing the function of a structural module, there are also the connection modules that describe how different structural modules are connected. Both are defined as `ModelingToolkit.jl` models. For this tutorial, we will limit ourselves to using functional modules already implemented in `PlantModules.jl`; implementation of new processes will be discussed in the second tutorial.
"""

# ╔═╡ 790b01e0-2736-4aec-8516-d22c04a22131
md"""
We list the pre-implemented functional modules here with a short description of each.
- Hydraulics
    - `hydraulic_module`: describes the hydraulics-driven growth of a compartment.
    - `environmental_module`: describes a compartment that is subject to hydraulic laws but does therefor not change in size.
- Carbon dynamics
    - `constant_carbon_module`: describes a compartment with a constant amount of carbon.
    - `simple_photosynthesis_module`: describes a compartment with carbon being produced during day and consumed proportionally to its concentration.
- Total water potential
    - `Ψ_soil_module`: describes an empirical relationship between the total water potential and relative water content of soil.
    - `Ψ_air_module`: describes the relationship between the total water potential and water potential of air.
- Hydraulic conductivity
    - `K_module`: describes hydraulic conductivity as proportional to the cross area of the compartment. 
    - `constant_K_module`: describes a constant hydraulic conductivity.
- Connection modules
    - `hydraulic_connection`: describes a water flow between two compartments based on the specific hydraulic conductivities of those compartments.
    - `daynight_hydraulic_connection`: describes a water flow between two compartments based on the hydraulic conductivities of those compartments that decreases during night.
    - `constant_hydraulic_connection`: describes a water flow between two compartments and specifies its own hydraulic conductivity.
"""

# ╔═╡ 1e2eaa86-d6e0-4749-bc92-04c64fe9f47d
md"""
Please note that the carbon dynamics modules are very simplistic and only exist so users can make a model that works out of the box. Users that want to create a realistic model are expected to implement some functional modules of their own in accordance with their own use for the package, as will be discussed next tutorial.
"""

# ╔═╡ 930e7ed8-0bfe-4e5a-8890-a1d1ce155881
md"""
#### Mapping function to structure
"""

# ╔═╡ 4cedbd9d-84ed-46f3-9a10-6cb993643f87
md"""
When the functional modules have been defined, we can assign them to the different types of structural module and the connections between them. We can do this for both with a simple `Dictionary`. Note that, for now, while structural modules can be assigned multiple functional modules, their connections can only be assigned one connection module.
"""

# ╔═╡ d54705b3-d8f4-4cc2-a780-369343749113
module_coupling = Dict(
	:Stem => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, simple_photosynthesis_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module]
);

# ╔═╡ 611289e9-e22c-4e6e-beec-ccea90eb48c9
connecting_modules = Dict(
	(:Soil, :Stem) => hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Leaf) => hydraulic_connection,
	(:Leaf, :Air) => daynight_hydraulic_connection,
);

# ╔═╡ b8709a33-301c-41e6-bb24-1035c68a7612
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ 4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
md"""
#### Defining parameter and initial values 
"""

# ╔═╡ 3035b6d0-bca0-4803-b32a-da1459bdd880
md"""
Next to the equations themselves, an important part in describing the plant's functional processes is defining the equations' parameter - and initial values. Since we'll restrict ourselves to hydraulic equations for this tutorial, these parameters and variables are the ones described in the section [Theoretical overview](nothinghere).

In PlantModules, parameter - and initial values are defined in a hierarchical manner, here listed in ascending priority:

> **Model-wide default values < module-specific default values < node-specific values**

For our current example, this means that there is no point in specifying the initial values for the compartment dimensions *D* for our Leaf compartments in the module-specific default values or the model-wide default values, since we already defined these values in all of the Leaf nodes of the graph, which have the highest priority.
"""

# ╔═╡ df4cd3de-d2b2-4f11-b755-a36e640fd2d5
md"""
Before we change any of the default parameter - and initial values, we can take a look at the defaults used by PlantModules:
"""

# ╔═╡ 1394e5ed-39f3-4a9f-8737-be7183219314
PlantModules.default_values

# ╔═╡ fbe00b62-6ca7-4c90-9050-081312911c74
md"""
Now imagine we have some information on our plant in question and it calls for different default values. These changes will be inputted as a Dictionary:
"""

# ╔═╡ 271d48a7-7022-4766-83d9-a70fab92515e
default_changes = Dict(:Ψ => -0.2);

# ╔═╡ 3c13c600-5135-44ea-8fc2-a1e11f72d0c5
md"""
Next to changing functional values over the entire model, some of the information we have may require different default values specifically for certain structural modules. We can specify this as shown below.
"""

# ╔═╡ 5f21a4b0-f663-4777-94f3-d00acba109b3
module_defaults = Dict(
	:Stem => Dict(:D => [0.5, 10.0], :M => 200e-6),
	:Leaf => Dict(:shape => Cuboid(), :M => 300e-6),
	:Soil => Dict(:W_max => 1e3, :K => 1.0),
	:Air => Dict(:W_r => 0.7, :K => 1e-3)
);

# ╔═╡ e0e5e7f7-d2f5-404b-a6c6-6848c318eccb
md"""
As you can see, we changed the default $shape$ for the leaf modules, as well as the initial dimensions $D$ and metabolite concentration $M$ of the stem. For both environmental compartments, we changed the hydraulic conductivity, and additionally the maximum water content for the soil comparment. For a detailed explanation of the parameters and initial values here, we again refer to the [theoretical overview](nothinghere).
"""

# ╔═╡ 53b79b18-524a-40b4-98e8-5900d7124305
md"""
Values for the connection modules can be modified in a similar way:
"""

# ╔═╡ 7d8a8a86-e40c-46fa-8cdd-8e15524f305f
connection_values = Dict(
	(:Leaf, :Air) => Dict(:t_sunrise => 6, :t_sunset => 22),
);

# ╔═╡ fb72735b-3d45-438d-ad83-3e36f42f5bb8
md"""
!!! note
	In contrast to the nodes, PlantModules currently does not support edges of the same type (that is, between the same type of structural module) to specify different parameter - or initial values in the graph definition.
"""

# ╔═╡ 113c0bdc-53e4-4a19-a47c-f4afba237eeb
md"""
Finally, all parameter changes are wrapped in the `PlantParameters` struct.
"""

# ╔═╡ ebf46963-bb9d-4604-924c-2b051189debc
plantparams = PlantParameters(; default_changes, module_defaults, connection_values);

# ╔═╡ 210d81ef-153e-4744-8266-80af4099770c
md"""
### Bringing it all together
"""

# ╔═╡ bc7573e7-bcd6-4347-9b0c-9111a824c9b5
md"""
Now that we have all parts of our model defined, all that's left to do is putting it together. For this we use the main workhorse of PlantModules: `generate_system`.
"""

# ╔═╡ a3c5dba8-8acc-424a-87bc-d1c6a263578c
system = generate_system(plantstructure, plantcoupling, plantparams);

# ╔═╡ d51795b2-32d3-455c-b785-5e501cfbdd08
md"""
This function will generate the `ODESystem` describing the model. It is possible to fine-tune the model even further at this stage as described in the [Customizing the model](nothinghere) section of the docs, thought this should generally not be required.
"""

# ╔═╡ d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
md"""
## Running the model
"""

# ╔═╡ 6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
md"""
Running the model is taken care of DifferentialEquations.jl. Users that are unfamiliar with the package may want to take a brief look at [the DifferentialEquations.jl docs](https://docs.sciml.ai/DiffEqDocs/stable/).
"""

# ╔═╡ bf114636-1e35-49f1-9407-f472b443a9ea
time_span = (0.0, 48.0);

# ╔═╡ 50d6fc31-80f5-4db7-b716-b26765008a0d
prob = ODEProblem(system, [], time_span, sparse = true);

# ╔═╡ c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
sol = solve(prob);

# ╔═╡ a399ea81-5392-4a54-8a40-953faf5b234a
md"""
In order to compare the results of higher transpiration, we remake and solve the problem with lower relative air humidity. This can be done by repeating the parameter, problem, system and solution definition with a different value for `W_r`, or we can simply use `ModelingToolkit.jl`'s `remake` function. To make it easier to get the correct variable from our system with many subsystems, the `get_subsystem_variables` convenience function exists that extract all variables corresponding to a certain name and structural module (or connection between two).
"""

# ╔═╡ 3526f17e-0f8b-4c10-9f33-96832673136d
air_W_r = get_subsystem_variables(system, plantstructure, :W_r, :Air)[1]

# ╔═╡ 63ce4478-3233-4216-9fa2-36fcf52e8673
prob2 = remake(prob, u0 = [air_W_r => 0.5]);

# ╔═╡ cb530432-97c8-4c1b-b6ef-905c1e1b5c81
sol2 = solve(prob2);

# ╔═╡ a6608eff-9399-443c-a33a-c62341f7b14c
md"""
## Showing the results
"""

# ╔═╡ a84c1948-26ac-497d-b2f5-8310e5341b52
md"First we verify that there is indeed a higher transpiration in our second simulation by plotting the total incoming water influx of the air."

# ╔═╡ f66ca207-98a2-40ee-bf95-ab6e191cc60f
plot(
	plotgraph(sol, plantstructure, varname = :ΣF, structmod = :Air, title = "Low transpiration", xaxis = "t (hr)"),
	plotgraph(sol2, plantstructure, varname = :ΣF, structmod = :Air, title = "High transpiration", xaxis = "t (hr)"), ylims = (0.0, 0.25), yaxis = "ΣF"
)

# ╔═╡ b523a45d-21b4-4bc1-9e47-70ebdb0c45f5
md"Finally, we can inspect whether the effect on the water potential of the stem is as expected. We plot the water potential of the leaves along with it to see if they follow the same pattern."

# ╔═╡ cf30d4f4-a5de-4def-8674-48088eabf17b
plot(
	plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Stem, :Leaf], title = "Low transpiration"),
	plotgraph(sol2, plantstructure, varname = :Ψ, structmod = [:Stem, :Leaf], title = "High transpiration"),
	ylims = (-0.35, -0.1), yaxis = "Ψ"
)

# ╔═╡ 79d012fd-4afd-4f3b-ad7c-8ca581bad1e5
md"""
The plot shows us the expected pattern for both plant parts, verifying that the pre-built functionality can capture this basic hydraulic phenomenon. For more interesting applications and more advanced functionality, we refer to the subsequent tutorials.
"""

# ╔═╡ a8a57cdb-65e7-4b9d-b8a5-f0b0ebf859dd
md"""
!!! note
	The water potential behaving differently on the second day compared to the first is caused by changes in the metabolite concentration between those days, illustrated in the plots below.
"""

# ╔═╡ 11212e1a-9bff-4ae2-8b88-b1bc0529543a
plotgraph(sol, plantstructure, varname = :M, structmod = [:Leaf, :Stem])

# ╔═╡ 1c6fac4d-7c63-4601-8b92-21196f5bb96b
plot(
	plotgraph(sol, plantstructure, varname = :Π, structmod = :Leaf),
	plotgraph(sol, plantstructure, varname = :P, structmod = :Leaf),
	plotgraph(sol, plantstructure, varname = :Ψ, structmod = :Leaf,
			  xlabel = "Time (hr)"),
	layout = (3, 1)
)

# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
# ╟─06004b15-17d4-4f43-a6f4-b67f5610ebac
# ╟─bc44fae7-6c4a-4f27-be51-dbabf4037601
# ╠═57b8dcb8-9baa-4ddf-9368-431b1be5850e
# ╠═662c0476-70aa-4a60-a81c-d7db2248b728
# ╠═65f88593-1180-447a-900f-49aef4647cd1
# ╠═7ff01837-8d59-4df7-ad50-08157224e166
# ╟─aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
# ╟─6ef5c63a-b753-43ae-baee-f6c24313a385
# ╟─1144887a-a4c7-46f6-9cf8-cba50cf873d0
# ╟─d131a541-9357-4d1a-b6a6-539bf5c9d884
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╟─aec7bcd6-6f27-4cf5-a955-f4d59e778fd3
# ╟─659e911c-8af2-4a66-855a-e333c41120c1
# ╟─0cc02e82-4fe8-4f27-a2d2-eb4bfba6b291
# ╠═6b7ebc68-f4a1-4ed6-b12b-e4ac5ee9b00a
# ╠═d57718c2-c77d-42a8-924f-ebdfcd51d919
# ╟─c4dc4961-ba2b-4b23-b80e-7d4eb8d9a9f4
# ╠═9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
# ╠═27160915-5943-4048-b363-d123846b02b8
# ╟─98eac4c4-b39a-4e11-917a-90b03d7385d1
# ╠═e00c5135-1d66-4dec-8283-40ebe06a8038
# ╠═dac02191-b640-40f5-a7d6-e6b06b946c23
# ╠═3bf42137-1551-44d6-b7af-eab13a97b6ef
# ╠═db8fe96d-c8c2-47c0-9377-281ce577a920
# ╟─c36668fa-3046-4967-a616-841a048ea7f9
# ╠═8980f1eb-e461-4624-9cce-83a7cf822349
# ╠═61bf737a-2226-42dc-b93a-a8f578048268
# ╟─20049311-d6e6-41d3-a0d8-8bad88c174f9
# ╟─668e02ee-ac78-4b3d-983d-402ec04584ef
# ╠═caab574a-05c5-4c0d-9ae4-19fd514a6c6c
# ╠═579c41a8-e054-474a-985e-314c9a7657ce
# ╟─f03a61ce-a0ff-43ff-abdd-2342f76e8c93
# ╟─43211f69-6bfe-4fd1-b474-65d0601558de
# ╟─a6552c16-b013-43af-9483-639a79944749
# ╟─52bb8235-abea-45c5-929c-5824bc8681c4
# ╟─c04564c4-4fb5-47bf-bc14-77aaebdece15
# ╟─790b01e0-2736-4aec-8516-d22c04a22131
# ╟─1e2eaa86-d6e0-4749-bc92-04c64fe9f47d
# ╟─930e7ed8-0bfe-4e5a-8890-a1d1ce155881
# ╟─4cedbd9d-84ed-46f3-9a10-6cb993643f87
# ╠═d54705b3-d8f4-4cc2-a780-369343749113
# ╠═611289e9-e22c-4e6e-beec-ccea90eb48c9
# ╠═b8709a33-301c-41e6-bb24-1035c68a7612
# ╟─4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
# ╟─3035b6d0-bca0-4803-b32a-da1459bdd880
# ╟─df4cd3de-d2b2-4f11-b755-a36e640fd2d5
# ╠═1394e5ed-39f3-4a9f-8737-be7183219314
# ╟─fbe00b62-6ca7-4c90-9050-081312911c74
# ╠═271d48a7-7022-4766-83d9-a70fab92515e
# ╟─3c13c600-5135-44ea-8fc2-a1e11f72d0c5
# ╠═5f21a4b0-f663-4777-94f3-d00acba109b3
# ╟─e0e5e7f7-d2f5-404b-a6c6-6848c318eccb
# ╟─53b79b18-524a-40b4-98e8-5900d7124305
# ╠═7d8a8a86-e40c-46fa-8cdd-8e15524f305f
# ╟─fb72735b-3d45-438d-ad83-3e36f42f5bb8
# ╟─113c0bdc-53e4-4a19-a47c-f4afba237eeb
# ╠═ebf46963-bb9d-4604-924c-2b051189debc
# ╟─210d81ef-153e-4744-8266-80af4099770c
# ╟─bc7573e7-bcd6-4347-9b0c-9111a824c9b5
# ╠═a3c5dba8-8acc-424a-87bc-d1c6a263578c
# ╟─d51795b2-32d3-455c-b785-5e501cfbdd08
# ╟─d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
# ╟─6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
# ╠═bf114636-1e35-49f1-9407-f472b443a9ea
# ╠═50d6fc31-80f5-4db7-b716-b26765008a0d
# ╠═c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
# ╟─a399ea81-5392-4a54-8a40-953faf5b234a
# ╠═3526f17e-0f8b-4c10-9f33-96832673136d
# ╠═63ce4478-3233-4216-9fa2-36fcf52e8673
# ╠═cb530432-97c8-4c1b-b6ef-905c1e1b5c81
# ╟─a6608eff-9399-443c-a33a-c62341f7b14c
# ╟─a84c1948-26ac-497d-b2f5-8310e5341b52
# ╠═f66ca207-98a2-40ee-bf95-ab6e191cc60f
# ╟─b523a45d-21b4-4bc1-9e47-70ebdb0c45f5
# ╠═cf30d4f4-a5de-4def-8674-48088eabf17b
# ╟─79d012fd-4afd-4f3b-ad7c-8ca581bad1e5
# ╟─a8a57cdb-65e7-4b9d-b8a5-f0b0ebf859dd
# ╠═11212e1a-9bff-4ae2-8b88-b1bc0529543a
# ╠═1c6fac4d-7c63-4601-8b92-21196f5bb96b
