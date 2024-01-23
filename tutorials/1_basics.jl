### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 34fd1713-4d0a-4bc9-81e1-bacf418747a2
using Pkg; Pkg.activate("..")

# ╔═╡ 65f88593-1180-447a-900f-49aef4647cd1
using ModelingToolkit, DifferentialEquations, Plots, Unitful, PlantGraphs # Should change to just this package and maybe VPL.jl stuff

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 1: Package basics
"""

# ╔═╡ 1144887a-a4c7-46f6-9cf8-cba50cf873d0
md"""
## Toy example
Consider a young potted plant growing on a windowsill indoors. We would like to simulate the evapotranspiration of water in the soil so we can estimate when the plant will need to watered again.
"""

# ╔═╡ 7d69ce40-8db3-4f11-bfae-695ececc94dc
md"""
### Loading the package.
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
myplant = Root() + Stem() + (Leaf([0.25, 0.65]), Leaf([0.15, 0.45]))

# ╔═╡ 43211f69-6bfe-4fd1-b474-65d0601558de
md"""
### 
"""

# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─1144887a-a4c7-46f6-9cf8-cba50cf873d0
# ╟─7d69ce40-8db3-4f11-bfae-695ececc94dc
# ╠═34fd1713-4d0a-4bc9-81e1-bacf418747a2
# ╠═65f88593-1180-447a-900f-49aef4647cd1
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╠═659e911c-8af2-4a66-855a-e333c41120c1
# ╠═e920f6aa-4c7b-4fd1-9dca-d9e3d4155ec2
# ╠═6b7ebc68-f4a1-4ed6-b12b-e4ac5ee9b00a
# ╠═d57718c2-c77d-42a8-924f-ebdfcd51d919
# ╟─a740d4ab-5ad8-4db4-9a80-aef2625a7d7b
# ╠═9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
# ╠═43211f69-6bfe-4fd1-b474-65d0601558de
