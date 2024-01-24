### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 34fd1713-4d0a-4bc9-81e1-bacf418747a2
using Pkg; Pkg.activate("..")

# ╔═╡ 65f88593-1180-447a-900f-49aef4647cd1
using PlantGraphs, GLMakie, ModelingToolkit #! MTK imports etc. should not be necessary when package is done

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 1: Package basics
"""

# ╔═╡ 1144887a-a4c7-46f6-9cf8-cba50cf873d0
md"""
## Toy example description
Consider a young potted plant growing on a windowsill indoors. We would like to simulate the evapotranspiration of water in the soil so we can estimate when the plant will need to watered again.
"""

# ╔═╡ 45481dd7-5c9c-4865-add9-9774d510c678
md"""
## Tackling the problem
"""

# ╔═╡ 7d69ce40-8db3-4f11-bfae-695ececc94dc
md"""
### Loading the package
"""

# ╔═╡ b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
md"""
### Defining the plant structure

The first step in modeling our plant's growth is defining its structure. There are a lot of options for formalizing the structure of a plant, and perhaps the most obvious first choice is what spatial scale to model on. PlantModules allows for a lot of freedom here, including the option to combine multiple spatial scales. 
For the first tutorial, we'll consider a very simple example: an organ-scale plant model of a plant with two leaves, a stem and a root system.
"""

# ╔═╡ 659e911c-8af2-4a66-855a-e333c41120c1
md"""
Regardless of the chosen structure, PlantModules expects it to be defined as a graph. We'll be using the graph implementation from the PlantGraphs.jl package for this tutorial, though any graph implementation can be used.
See [the PlantGraphs](https://virtualplantlab.com/stable/manual/Graphs/) docs for more information on the use of these graphs, and [the PlantModules](nothinghere) docs for more information on using custom graph implementations.
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
Individual graph nodes can contain parameter values and initial values to allow differences between nodes of the same type. Here, we'll give the leafs a size field so we can start off one of them larger than the other.
"""

# ╔═╡ 9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
plantstruc = Root() + Stem() + (Leaf([0.25]), Leaf([0.15]))

# ╔═╡ ecb23f1e-ee39-4c5a-911f-eaa252c23968
md"""
We can inspect the structure of our plant by visualising the graph:
"""

# ╔═╡ 86d81fa3-bea4-40fa-9522-7db9fe2f6a82
draw(plantstruc, resolution = (500, 400))

# ╔═╡ 98eac4c4-b39a-4e11-917a-90b03d7385d1
md"""
### Defining the environment

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
### Coupling function to structure

Now that we know what our plant is structured like, we need to tell PlantModules how those structural components behave functionally. For example, some plant parts may be photosynthetically active while others are not. In this tutorial, we'll mostly neglect carbon dynamics and instead focus on the plant's water dynamics.
"""

# ╔═╡ 3035b6d0-bca0-4803-b32a-da1459bdd880
md"""
Defining the plant's functional properties comes down to specifying values for the biophysical parameters of the equations describing the plant's growth, as well as initial values for the equations' variables. Since we'll restrict ourselves to hydraulic equations for this tutorial, these parameters and variables are the ones described in the section [Theoretical overview](nothinghere).

In PlantModules, parameter - and initial values are defined in a hierarchical manner, here listed in ascending priority:

> **Model-wide default values < module-specific default values < node-specific values**

For our current example, this means that there is no point in specifying the initial values for the compartment dimensions *D* for our Leaf compartments in the module-specific default values, since we already defined these values in the Leaf nodes of the graph.
"""

# ╔═╡ b282eb11-05c6-4985-aa66-2174a82ebb6e
# (from PlantModules, not written out by user)
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

# ╔═╡ c04564c4-4fb5-47bf-bc14-77aaebdece15
md"""
When coupling function to structure, PlantModules (as the name suggests) makes use of **modules** to describe the functional behaviour of a certain compartment type. These modules are sets of equations implemented in [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/). The user is free to implement these themselves, though here we will use the two already provided by PlantModules:
- hydraulic_module, which describes the hydraulics-driven growth of a compartment.
- environmental_module, which describes a compartment that is subject to hydraulic laws but does therefor not change in size.

The latter module, which we will use to simulate the soil and air, still expects the user to provide some equation describing how the water potential changes. We will define these first for both air and the soil.
"""

# ╔═╡ 68fbfd88-b1a6-4d52-aee4-37e76b191fe4
Ψ_soil_func(W_r) = -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19) # An empirical relationship between the soil water potential and relative water content

# ╔═╡ 10ac5d18-8527-4bb8-aa5d-0da752c9a808
md"""
The function for the air water potential uses two factors not yet described in the system: the ideal gas constant *R* and the molar volume of water *V_w*. We will define these in the global scope using the standard ModelingToolkit.jl syntax. 
"""

# ╔═╡ 137a367e-7b0e-4ef2-8068-628158f3a45d
@constants R = 8.314

# ╔═╡ 3a876d6a-0ea6-44ad-89be-7eac04acb7bb
@parameters V_w = 18e-6

# ╔═╡ b69ee1cb-6506-4152-9ef0-b02a43a90990
Ψ_air_func(T, W_r) = R * T / V_w * log(W_r)

# ╔═╡ a3c5dba8-8acc-424a-87bc-d1c6a263578c
plantsys = PlantSystem(
	structure = plantstruc,
	modules = [
		:hydraulic_module => [:Root, :Stem, :Leaf],
		:environmental_module => [:Soil, :Air]
	],
	module_info = [
		:Root => [:D => [0.3, 0.05, 0.03], :ϵ_D => [5.0, 0.3, 0.2], :ϕ_D => [0.7, 0.1, 0.05]],
		:Stem => [:D => [0.4, 0.03], :ϵ_D => [6.0, 0.15], :ϕ_D => [0.8, 0.03]],
		:Leaf => [:ϵ_D => [3.0], :ϕ_D => [0.45]],
		:Soil => [:W_max => 500.0, :T => 288.15, :Ψ => Ψ_soil_func],
		:Air => [:Ψ => Ψ_air_func]
	],
	connection_info = [
		(:Root, :Stem) => [:K => 800]
	],
	default_values = [:Γ => 0.4, :P => 0.1, :M => 15.0, :T => 293.15]
)

# ╔═╡ fb3c58df-1d6b-4ced-803d-2d0fc537b942
md"""
# TO-DO list

- Update hyperlinks with pages of docs when they exist
- Change parameter - and initial values with logical ones based on some data or something
- Add shapes in some better way
- Add Unitful units?
- How to add multiple soil compartments?
"""

# ╔═╡ 38c69eea-a4dd-4fc0-951f-dc36e9530b80


# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─1144887a-a4c7-46f6-9cf8-cba50cf873d0
# ╟─45481dd7-5c9c-4865-add9-9774d510c678
# ╟─7d69ce40-8db3-4f11-bfae-695ececc94dc
# ╠═34fd1713-4d0a-4bc9-81e1-bacf418747a2
# ╠═65f88593-1180-447a-900f-49aef4647cd1
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╟─659e911c-8af2-4a66-855a-e333c41120c1
# ╠═e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
# ╠═6b7ebc68-f4a1-4ed6-b12b-e4ac5ee9b00a
# ╠═d57718c2-c77d-42a8-924f-ebdfcd51d919
# ╟─a740d4ab-5ad8-4db4-9a80-aef2625a7d7b
# ╠═9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
# ╟─ecb23f1e-ee39-4c5a-911f-eaa252c23968
# ╠═86d81fa3-bea4-40fa-9522-7db9fe2f6a82
# ╟─98eac4c4-b39a-4e11-917a-90b03d7385d1
# ╠═e00c5135-1d66-4dec-8283-40ebe06a8038
# ╠═dac02191-b640-40f5-a7d6-e6b06b946c23
# ╟─43211f69-6bfe-4fd1-b474-65d0601558de
# ╟─3035b6d0-bca0-4803-b32a-da1459bdd880
# ╟─b282eb11-05c6-4985-aa66-2174a82ebb6e
# ╟─c04564c4-4fb5-47bf-bc14-77aaebdece15
# ╠═68fbfd88-b1a6-4d52-aee4-37e76b191fe4
# ╠═b69ee1cb-6506-4152-9ef0-b02a43a90990
# ╟─10ac5d18-8527-4bb8-aa5d-0da752c9a808
# ╠═137a367e-7b0e-4ef2-8068-628158f3a45d
# ╠═3a876d6a-0ea6-44ad-89be-7eac04acb7bb
# ╠═a3c5dba8-8acc-424a-87bc-d1c6a263578c
# ╠═fb3c58df-1d6b-4ced-803d-2d0fc537b942
# ╠═38c69eea-a4dd-4fc0-951f-dc36e9530b80
