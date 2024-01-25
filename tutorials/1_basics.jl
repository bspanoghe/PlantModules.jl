### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 34fd1713-4d0a-4bc9-81e1-bacf418747a2
# Maybe not include this in the tutorial?
using Pkg; Pkg.activate("..")

# ╔═╡ 65f88593-1180-447a-900f-49aef4647cd1
using PlantGraphs, GLMakie, ModelingToolkit #! MTK imports etc. should not be necessary when package is done

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 1: Package basics
"""

# ╔═╡ d16a6d55-1f29-4b98-b1af-2dee1d38f386
module PlantModules

export PlantSystem

function hydraulic_module(; name) end

function environmental_module(; name) end

struct PlantSystem
	structure
	modules::Vector{Pair{Symbol, Vector{Symbol}}}
	module_info::Vector{Pair{Symbol, Vector}}
	connection_info::Vector
	default_values::Vector{Pair{Symbol, T}} where {T}

	function PlantSystem(; structure, modules, module_info, connection_info, default_values)
		#! calculate MTK things here and add to struct
		return new(structure, modules, module_info, connection_info, default_values)
	end
end

end # module

# ╔═╡ 6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
md"""
## Introduction 👋

In this tutorial, we'll cover the most basic functionalities of the PlantModules package, which is the following:
- 🐍
"""

# ╔═╡ 1144887a-a4c7-46f6-9cf8-cba50cf873d0
md"""
### Toy example description

For the first tutorial, let's start with a very simple problem. Consider a young potted plant growing on a windowsill indoors. We would like to simulate the evapotranspiration of water in the soil so we can estimate when the plant will need to be watered again.

![plantfigu](https://www.almanac.com/sites/default/files/users/The%20Editors/pepper-seedlings_full_width.jpg)
"""
#! replace picture with something royalty free. a drawing?

# ╔═╡ aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
md"""
## Creating the model 🛠
"""

# ╔═╡ 6ef5c63a-b753-43ae-baee-f6c24313a385
md"""
As the name suggests, one of the most important goals of PlantModules is to enable the user to easily model plant growth in a _modular_ manner. In order to achieve this, we will define our plant in function of a few sets of similarly behaving parts or "modules" for short. These modules can be similar both on the structural - and on the functional level.

On the structural level, some examples are:
- an oak tree can be considered a repeating module in a forest.
- a branch can be considered a repeating module in a tree.
- a collenchyma cell can be considered a repeating module in a branch.

On the functional level, on the other hand:
- most structural modules of a tree (branches, leaves, root segments, etc.) contain water which will flow between them as dictated by hydraulic laws, which plays an important role in plant growth.
- one or more of the plant's structural modules (notably the leaves) will assimilate carbon through photosynthesis, while the others do not.

As such, the workflow to create a model in PlantModules boils down to defining these modules and how they interact. Let's jump right in!
"""

# ╔═╡ b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
md"""
### Defining the structural modules
"""

# ╔═╡ aec7bcd6-6f27-4cf5-a955-f4d59e778fd3
md"""
#### The plant

The first step in modeling our plant's growth is defining its structure. There are a lot of options for formalizing the structure of a plant, and perhaps the most obvious first choice is what spatial scale to model on. PlantModules allows for a lot of freedom here, including the option to combine multiple spatial scales. 
For the first tutorial, we'll consider a very simple example: an organ-scale plant model of a plant with two leaves, a stem and a root system.
"""

# ╔═╡ 659e911c-8af2-4a66-855a-e333c41120c1
md"""
Regardless of the chosen structure, PlantModules expects it to be defined as a graph. We'll be using the graph implementation from the PlantGraphs.jl package for this tutorial, though any graph implementation can be used.
See [the PlantGraphs](https://virtualplantlab.com/stable/manual/Graphs/) docs for more information on the use of these graphs, and [the PlantModules](nothinghere) docs for more information on using custom graph implementations.
"""

# ╔═╡ e232199f-ee2f-4294-8762-f41b37883d26
md"""
As any good package tutorial, we start by loading in the required packages.
"""

# ╔═╡ 0cc02e82-4fe8-4f27-a2d2-eb4bfba6b291
md"""
Then we can define the functional modules of our plant. For our example, there are only three.
"""

# ╔═╡ e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
mutable struct Root <: Node end

# ╔═╡ 6b7ebc68-f4a1-4ed6-b12b-e4ac5ee9b00a
mutable struct Stem <: Node end

# ╔═╡ d57718c2-c77d-42a8-924f-ebdfcd51d919
mutable struct Leaf <: Node
	D::Vector
end

# ╔═╡ a740d4ab-5ad8-4db4-9a80-aef2625a7d7b
md"""
Individual graph nodes can contain parameter values and initial values to allow differences between nodes of the same type. Here, we'll give the leafs a size field so we can start off one of them larger than the other. We'll see more on parameter values later.
"""

# ╔═╡ b0b5e539-19f5-4feb-8a2f-a6b2f0590c7c
md"""
To continue, we specify how our structural modules interact with each other.
"""

# ╔═╡ 9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
plantstruc = Root() + Stem() + (Leaf([0.25]), Leaf([0.15]))

# ╔═╡ ecb23f1e-ee39-4c5a-911f-eaa252c23968
md"""
We can inspect the structure of our plant by visualising the graph:
"""

# ╔═╡ 86d81fa3-bea4-40fa-9522-7db9fe2f6a82
draw(plantstruc, resolution = (500, 400))

# ╔═╡ f03a61ce-a0ff-43ff-abdd-2342f76e8c93
md"""
For our tiny example plant, we only have one structural module that actually repeats, somewhat defeating the purpose of modelling in a modular manner: we may as well write this entire model out by hand! Rest assured, however, that the approach we're seeing here also works for trees with hundreds of branches and thousands of leaves, as we will see in the next tutorial.
"""

# ╔═╡ 98eac4c4-b39a-4e11-917a-90b03d7385d1
md"""
#### The environment

Plants need an environment to grow in. For most plants, the most basic environmental compartments that need defining are **the soil** from which the plant gets water and nutrients, and **the air** with which the plant exchanges gasses and to which it loses water. 

For more intricate models one may, for example, also want to take into account the sun as well as possible sources of shadow or divide the soil into multiple compartments, but we will stick to the basics for now.
"""

# ╔═╡ e00c5135-1d66-4dec-8283-40ebe06a8038
struct Soil <: Node
	W
end

# ╔═╡ dac02191-b640-40f5-a7d6-e6b06b946c23
struct Air <: Node end

# ╔═╡ 43211f69-6bfe-4fd1-b474-65d0601558de
md"""
### Defining the functional modules

Now that we know what our plant is structured like, we need to tell PlantModules how those structural components behave functionally. Some examples, such as only a select set of structural modules doing photosynthesis, have been mentioned before. In this tutorial, we'll mostly neglect carbon dynamics and instead focus on the plant's water dynamics.
"""

# ╔═╡ c04564c4-4fb5-47bf-bc14-77aaebdece15
md"""
PlantModules defines functional modules as as sets of differential equations implemented in [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/).
There are two functional modules already provided by PlantModules:

- `PlantModules.hydraulic_module`, which describes the hydraulics-driven growth of a compartment.
- `PlantModules.environmental_module`, which describes a compartment that is subject to hydraulic laws but does therefor not change in size.

The user is free to extend these or implement new ones themselves, though for now we will stick to using these two.
"""

# ╔═╡ 2ca2f739-2b61-4519-bdc4-d3081c793446
md"""
While these pre-implemented modules already define most of the functional process, they still expect some part to be defined by the user.
The former module, which we will use to simulate the three structural plant modules, requires an equation defining the carbon content over time.
The latter module, which we will use to simulate the soil and air, still expects the some equation describing how the water potential changes. 

We will define these equations first. Note that you can use any of the variables and parameters defined in the functional module.
"""
#! Replace with default behaviour and skip here?

# ╔═╡ 96064ce5-d555-46a7-a647-8f94de01cd31
C_root_func() = 0.5 # We'll assume the soluble carbon content remains constant over the simulated time 

# ╔═╡ 805e3bf4-0b17-4a7e-a9e6-b18d64ea52bb
C_stem_func() = 1

# ╔═╡ af203a2b-2299-4b2b-b034-c3cb39648bb7
C_leaf_func() = 3

# ╔═╡ 68fbfd88-b1a6-4d52-aee4-37e76b191fe4
Ψ_soil_func(W_r) = -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19) # An empirical relationship between the soil water potential and relative water content

# ╔═╡ 10ac5d18-8527-4bb8-aa5d-0da752c9a808
md"""
The function for the air water potential uses two constants not yet described in the system: the ideal gas constant *R* and the molar volume of water *V_w*. We will define these in the global scope using the standard ModelingToolkit.jl syntax. 
"""

# ╔═╡ 137a367e-7b0e-4ef2-8068-628158f3a45d
@constants R = 8.314 V_w = 18e-6

# ╔═╡ b69ee1cb-6506-4152-9ef0-b02a43a90990
Ψ_air_func(T, W_r) = R * T / V_w * log(W_r)

# ╔═╡ 930e7ed8-0bfe-4e5a-8890-a1d1ce155881
md"""
### Coupling functional and structural modules
"""

# ╔═╡ 4cedbd9d-84ed-46f3-9a10-6cb993643f87
md"""
Our model needs to know which structural modules make use of which functional modules. As perhaps the simplest part of our modeling workflow, we can define this using a simple vector of pairs.
"""

# ╔═╡ d54705b3-d8f4-4cc2-a780-369343749113
module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air]
]

# ╔═╡ 4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
md"""
### Defining parameter and initial values 
"""

# ╔═╡ 3035b6d0-bca0-4803-b32a-da1459bdd880
md"""
Next to the equations themselves, an important part in describing the plant's functional processes is defining the equations' parameter - and initial values. Since we'll restrict ourselves to hydraulic equations for this tutorial, these parameters and variables are the ones described in the section [Theoretical overview](nothinghere).

In PlantModules, parameter - and initial values are defined in a hierarchical manner, here listed in ascending priority:

> **Model-wide default values < module-specific default values < node-specific values**

For our current example, this means that there is no point in specifying the initial values for the compartment dimensions *D* for our Leaf compartments in the module-specific default values, since we already defined these values in the Leaf nodes of the graph, which have the highest priority.
"""

# ╔═╡ 5f21a4b0-f663-4777-94f3-d00acba109b3
module_defaults = [
	:Root => [:D => [0.3, 0.05, 0.03], :ϵ_D => [5.0, 0.3, 0.2], :ϕ_D => [0.7, 0.1, 0.05]],
	:Stem => [:D => [0.4, 0.03], :ϵ_D => [6.0, 0.15], :ϕ_D => [0.8, 0.03]],
	:Leaf => [:ϵ_D => [3.0], :ϕ_D => [0.45]],
	:Soil => [:W_max => 500.0, :T => 288.15, :Ψ => Ψ_soil_func],
	:Air => [:Ψ => Ψ_air_func]
]

# ╔═╡ 271d48a7-7022-4766-83d9-a70fab92515e
model_defaults = [:Γ => 0.4, :P => 0.1, :M => 15.0, :T => 293.15]

# ╔═╡ 611289e9-e22c-4e6e-beec-ccea90eb48c9
myconnection_info = [
	(:Root, :Stem) => [:K => 800],
	(:Leaf, :Air) => [:K => 1e-3]
]

# ╔═╡ 210d81ef-153e-4744-8266-80af4099770c
md"""
### Bringing it all together
"""

# ╔═╡ a3c5dba8-8acc-424a-87bc-d1c6a263578c
plantsys = PlantModules.PlantSystem(
	structure = plantstruc,
	modules = module_coupling,
	module_info = module_defaults,
	connection_info = myconnection_info,
	default_values = model_defaults
)

# ╔═╡ fb3c58df-1d6b-4ced-803d-2d0fc537b942
md"""
# TO-DO list

- Update hyperlinks with pages of docs when they exist
- Change parameter - and initial values with logical ones based on some data or something
- Add shapes in some better way
- Add Unitful units?
- How to add multiple soil compartments?
- Make terminology for plant part types / compartments / modules more consistent
- Explain connections better! Also allow for custom connections?
"""

# ╔═╡ 38c69eea-a4dd-4fc0-951f-dc36e9530b80


# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─d16a6d55-1f29-4b98-b1af-2dee1d38f386
# ╟─6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
# ╟─1144887a-a4c7-46f6-9cf8-cba50cf873d0
# ╟─aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
# ╟─6ef5c63a-b753-43ae-baee-f6c24313a385
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╟─aec7bcd6-6f27-4cf5-a955-f4d59e778fd3
# ╟─659e911c-8af2-4a66-855a-e333c41120c1
# ╟─e232199f-ee2f-4294-8762-f41b37883d26
# ╟─34fd1713-4d0a-4bc9-81e1-bacf418747a2
# ╠═65f88593-1180-447a-900f-49aef4647cd1
# ╟─0cc02e82-4fe8-4f27-a2d2-eb4bfba6b291
# ╠═e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
# ╠═6b7ebc68-f4a1-4ed6-b12b-e4ac5ee9b00a
# ╠═d57718c2-c77d-42a8-924f-ebdfcd51d919
# ╟─a740d4ab-5ad8-4db4-9a80-aef2625a7d7b
# ╟─b0b5e539-19f5-4feb-8a2f-a6b2f0590c7c
# ╠═9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
# ╟─ecb23f1e-ee39-4c5a-911f-eaa252c23968
# ╠═86d81fa3-bea4-40fa-9522-7db9fe2f6a82
# ╟─f03a61ce-a0ff-43ff-abdd-2342f76e8c93
# ╟─98eac4c4-b39a-4e11-917a-90b03d7385d1
# ╠═e00c5135-1d66-4dec-8283-40ebe06a8038
# ╠═dac02191-b640-40f5-a7d6-e6b06b946c23
# ╟─43211f69-6bfe-4fd1-b474-65d0601558de
# ╟─c04564c4-4fb5-47bf-bc14-77aaebdece15
# ╟─2ca2f739-2b61-4519-bdc4-d3081c793446
# ╠═96064ce5-d555-46a7-a647-8f94de01cd31
# ╠═805e3bf4-0b17-4a7e-a9e6-b18d64ea52bb
# ╠═af203a2b-2299-4b2b-b034-c3cb39648bb7
# ╠═68fbfd88-b1a6-4d52-aee4-37e76b191fe4
# ╠═b69ee1cb-6506-4152-9ef0-b02a43a90990
# ╟─10ac5d18-8527-4bb8-aa5d-0da752c9a808
# ╠═137a367e-7b0e-4ef2-8068-628158f3a45d
# ╟─930e7ed8-0bfe-4e5a-8890-a1d1ce155881
# ╟─4cedbd9d-84ed-46f3-9a10-6cb993643f87
# ╠═d54705b3-d8f4-4cc2-a780-369343749113
# ╟─4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
# ╟─3035b6d0-bca0-4803-b32a-da1459bdd880
# ╠═5f21a4b0-f663-4777-94f3-d00acba109b3
# ╠═271d48a7-7022-4766-83d9-a70fab92515e
# ╠═611289e9-e22c-4e6e-beec-ccea90eb48c9
# ╟─210d81ef-153e-4744-8266-80af4099770c
# ╠═a3c5dba8-8acc-424a-87bc-d1c6a263578c
# ╟─fb3c58df-1d6b-4ced-803d-2d0fc537b942
# ╠═38c69eea-a4dd-4fc0-951f-dc36e9530b80
