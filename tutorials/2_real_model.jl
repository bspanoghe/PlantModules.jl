### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 34fd1713-4d0a-4bc9-81e1-bacf418747a2
# Maybe not include this in the tutorial?
using Pkg; Pkg.activate("..")

# ╔═╡ 65f88593-1180-447a-900f-49aef4647cd1
using PlantGraphs, GLMakie, ModelingToolkit, DifferentialEquations #! MTK imports etc. should not be necessary when package is done

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 2: Making a real model
"""

# ╔═╡ d16a6d55-1f29-4b98-b1af-2dee1d38f386
module PlantModules
#! To be actually written out so this thing works

export PlantSystem
import ModelingToolkit.ODEProblem, Plots.plot, DifferentialEquations.solve

function hydraulic_module(; name) end

function environmental_module(; name) end

function hydraulic_connection(; name) end

struct PlantSystem
	model_defaults
	module_defaults
	module_coupling
	struct_connections
	func_connections
	MTK

	function PlantSystem(; model_defaults, module_defaults, module_coupling, struct_connections, func_connections)
		MTK = missing #! calculate MTK things here and add to struct
		return new(model_defaults, module_defaults, module_coupling, struct_connections, func_connections, MTK)
	end
end

default_vals = [:T => 298.15, :P => 0.1, :Γ => 0.3]

function ODEProblem(ps::PlantSystem, timespan::Tuple) end

function solve(prob) end #! this one doesn't actually have to be extended

function plot(sol, struct_modules::Vector, func_vars::Vector)  end

end # module

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

There is no fun in modeling without purpose, which is why we'll introduce the package with a simple problem.

Our recently sprouted pepper seedlings are growing on the windowsill indoors. We just watered them, and we'd like to know how the water content of the soil changes through time so we can water them again at the optimal moment: the soil should have dried adequately to promote root growth without putting the seedlings under too much drought stress.

![plantfigu](https://www.almanac.com/sites/default/files/users/The%20Editors/pepper-seedlings_full_width.jpg)
"""
#! replace picture with something royalty free. a drawing?

# ╔═╡ aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
md"""
## Creating the model 🛠
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

As such, the workflow to create a model in PlantModules boils down to defining these modules and how they interact. Let's jump right in!
"""

# ╔═╡ b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
md"""
### Defining the structural modules
"""

# ╔═╡ aec7bcd6-6f27-4cf5-a955-f4d59e778fd3
md"""
#### The plant

The first step in modeling our plant's growth is defining its structure. There are a lot of options for formalizing the structure of a plant, and perhaps the most obvious first choice we need to make is what spatial scale to model on. PlantModules allows for a lot of freedom here, including the option to combine multiple spatial scales. 
For the first tutorial, we'll consider a very simple example: an organ-scale plant model of a plant with two leaves, a stem and a root system.
"""

# ╔═╡ 659e911c-8af2-4a66-855a-e333c41120c1
md"""
PlantModules uses **graphs** to define the relationships between modules. As such, structural modules need to be implemented as graph nodes of some sort. We'll be using the graph implementation from [the PlantGraphs](https://virtualplantlab.com/stable/manual/Graphs/) package for this tutorial, though [any graph implementation can be used](nothinghere).
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
Individual graph nodes can contain parameter values and initial values to allow differences between nodes of the same type. Here, we'll give the leaves a size field `D` so we can start off one of them larger than the other. We'll see more on parameter - and initial values later.
"""

# ╔═╡ 98eac4c4-b39a-4e11-917a-90b03d7385d1
md"""
#### The environment

Plants need an environment to grow in. For most plants, the most basic environmental compartments that need defining are the soil, from which the plant gets water and nutrients, and the air, with which the plant exchanges gasses and to which it loses water. 

For more intricate models one may for example also want to take into account the sun as well as possible sources of shadow or divide the soil into multiple compartments, but we will stick to the basics for now.
"""

# ╔═╡ e00c5135-1d66-4dec-8283-40ebe06a8038
struct Soil <: Node end

# ╔═╡ dac02191-b640-40f5-a7d6-e6b06b946c23
struct Air <: Node end

# ╔═╡ 43211f69-6bfe-4fd1-b474-65d0601558de
md"""
### Defining the functional modules

Now that we have an idea of the plant's structural modules, we need to define some sort of functionality expressed by one or more of them. Some examples, such as only a select set of structural modules doing photosynthesis, have been mentioned before. In this tutorial, we'll mostly neglect carbon dynamics and instead focus on the plant's water dynamics.
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
While these pre-implemented modules already define most of the functional process, they both lack an equation for one of the variables that they still expect to be defined by the user.
The former module, which we will use to simulate the three structural plant modules, requires a specification for the carbon content (in mol / m^3).
The latter module, which we will use to simulate the soil and air, still expects a specification of the water potential (in MPa).

We can provide this missing information as either a constant value or a function of the variables and parameters defined in that functional module. We'll define them here and connect them to their function modules right after.
"""
#! Replace with default behaviour and skip here?

# ╔═╡ 96064ce5-d555-46a7-a647-8f94de01cd31
C_root = 0.5 # We'll assume the soluble carbon content remains constant over the simulated time

# ╔═╡ 805e3bf4-0b17-4a7e-a9e6-b18d64ea52bb
C_stem = 1 # Idem here

# ╔═╡ af203a2b-2299-4b2b-b034-c3cb39648bb7
C_leaf = 3 # And here!

# ╔═╡ 68fbfd88-b1a6-4d52-aee4-37e76b191fe4
Ψ_soil_func(W_r) = -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19) # An empirical relationship between the soil water potential and relative water content

# ╔═╡ b69ee1cb-6506-4152-9ef0-b02a43a90990
Ψ_air_func(T, W_r) = R * T / V_w * log(W_r) #! What's this equation called again? (Ask Jeroen)

# ╔═╡ 10ac5d18-8527-4bb8-aa5d-0da752c9a808
md"""
The function for the air water potential uses two constants not yet described in the system: the ideal gas constant *R* and the molar volume of water *V_w*. We will define these in the global scope using the standard ModelingToolkit.jl syntax. 
"""

# ╔═╡ 137a367e-7b0e-4ef2-8068-628158f3a45d
@constants R = 8.314 V_w = 18e-6

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
PlantModules.default_vals

# ╔═╡ fbe00b62-6ca7-4c90-9050-081312911c74
md"""
Now imagine we have some data on our plant in question and it calls for different default values. This is done by specifying the name of the parameter or variable together with its desired value as a Julia Pair.
"""

# ╔═╡ 271d48a7-7022-4766-83d9-a70fab92515e
model_defaults = [:Γ => 0.4, :P => 0.2, :M => 15.0, :T => 293.15]

# ╔═╡ 3c13c600-5135-44ea-8fc2-a1e11f72d0c5
md"""
Some of the information we have requires different default values specifically for certain structural modules. We can specify this as shown below.
"""

# ╔═╡ 5f21a4b0-f663-4777-94f3-d00acba109b3
module_defaults = [
	:Root => [:D => [0.3, 0.05, 0.03], :ϵ_D => [5.0, 0.3, 0.2], :ϕ_D => [0.7, 0.1, 0.05], :C => C_root],
	:Stem => [:D => [0.4, 0.03], :ϵ_D => [6.0, 0.15], :ϕ_D => [0.8, 0.03], :C => C_stem],
	:Leaf => [:ϵ_D => [3.0], :ϕ_D => [0.45], :C => C_leaf],
	:Soil => [:W_max => 500.0, :T => 288.15, :Ψ => Ψ_soil_func],
	:Air => [:Ψ => Ψ_air_func]
]

# ╔═╡ d574346b-e8d8-49cd-a3bc-bf53617119f1
md"""
Notice that instead of a constant value, we specified a **function** for some of our variables, as discussed in the section [Defining the functional modules](nothinghere). The user is not limited to only providing functions for variables missing an equation in their functional module, though: functions can also be provided for parameters in case time-variable parameters are desired!
"""

# ╔═╡ 930e7ed8-0bfe-4e5a-8890-a1d1ce155881
md"""
### Coupling functional and structural modules
"""

# ╔═╡ 4cedbd9d-84ed-46f3-9a10-6cb993643f87
md"""
Our model needs to know which structural modules make use of which functional modules. As perhaps the simplest part of our modeling workflow, we can define this using another Vector of Pairs.
"""

# ╔═╡ d54705b3-d8f4-4cc2-a780-369343749113
module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air]
]

# ╔═╡ 5885b970-cfb6-4e3b-8039-23a5f7ab4bae
md"""
### Module connections
"""

# ╔═╡ b0b5e539-19f5-4feb-8a2f-a6b2f0590c7c
md"""
The final information required to fully define our plant growth model is how the modules relate to each other. This can once again be separated into a structural and functional part.
"""

# ╔═╡ c4d42004-20dd-40ec-9200-1dd2b69b0bec
md"""
#### Structural connections

As mentioned near the beginning, PlantModules uses graphs to define relationships between modules, and the structural modules were implemented as graph node types. All we need to do now is create a graph using these node types reflecting the desired structure of our plant.

Since we're using the PlantGraphs Julia package for this tutorial, we can easily define simple graphs using the DSL:
"""

# ╔═╡ 9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
plant_graph = Root() + Stem() + (Leaf([0.25]), Leaf([0.15]))

# ╔═╡ ecb23f1e-ee39-4c5a-911f-eaa252c23968
md"""
We can inspect the structure of our plant by visualising the graph:
"""

# ╔═╡ 86d81fa3-bea4-40fa-9522-7db9fe2f6a82
draw(plant_graph, resolution = (500, 400))

# ╔═╡ c36668fa-3046-4967-a616-841a048ea7f9
md"""
There are some structural modules missing here: the environmental ones. We could have added the soil compartment in the graph we just defined, but not the air. That is because both leaves are connected to this same air node, which would introduce a cycle into our graph, something many plant-related graph implementations do not support.

To overcome this, PlantModules allows the user to input multiple graphs and the ability to define what structural modules should be fully connected between all input graphs to acquire the final structural graph:
"""

# ╔═╡ 3bf42137-1551-44d6-b7af-eab13a97b6ef
soil_graph = Soil()

# ╔═╡ db8fe96d-c8c2-47c0-9377-281ce577a920
air_graph = Air()

# ╔═╡ f267ca25-5b7f-48e4-9865-1c7952215a78


# ╔═╡ 61bf737a-2226-42dc-b93a-a8f578048268
intergraph_connections = [(:Air, :Leaf), (:Soil, :Root)]

# ╔═╡ 668e02ee-ac78-4b3d-983d-402ec04584ef
md"""
The full structural connection information can then be acquired by putting them together in a vector.
"""

# ╔═╡ caab574a-05c5-4c0d-9ae4-19fd514a6c6c
struct_connections = [plant_graph, soil_graph, air_graph, intergraph_connections]

# ╔═╡ f03a61ce-a0ff-43ff-abdd-2342f76e8c93
md"""
> For our tiny example plant, we only have one structural module that actually repeats, somewhat defeating the purpose of modelling in a modular manner: we may as well write this entire model out by hand! Rest assured, however, that the approach we're seeing here also works for trees with hundreds of branches and thousands of leaves, as we will see in the next tutorial.
"""

# ╔═╡ ac6e098d-98ff-4940-8386-dc76c75eb2c3
md"""
#### Functional connections

If we consider functional modules as the information stored in the nodes of the graph, then the functional connections are the information stored in the edges. For the functional process under consideration here, this corresponds with a ModelingToolkit equation set describing water flow between compartments. It only has one parameter *K*, which is the hydraulic conductivity of the connection in question. 

Defining a new functional connection module is discussed in the following tutorial, for now we'll simply change the parameter value of this one:
"""

# ╔═╡ 611289e9-e22c-4e6e-beec-ccea90eb48c9
func_connections = [
	:default => PlantModules.hydraulic_connection,
	(:Soil, :Root) => [:K => 50],
	(:Root, :Stem) => [:K => 800],
	(:Leaf, :Air) => [:K => 1e-3]
]

# ╔═╡ fb72735b-3d45-438d-ad83-3e36f42f5bb8
md"""
> ⚠ In contrast to the nodes, PlantModules currently does not support edges of the same type (that is, between the same structural modules) to specify different parameter - or initial values in the graph definition.
"""

# ╔═╡ 210d81ef-153e-4744-8266-80af4099770c
md"""
### Bringing it all together
"""

# ╔═╡ bc7573e7-bcd6-4347-9b0c-9111a824c9b5
md"""
Now that we have all parts of our model defined, all that's left to do is putting it together. For this we use the main workhorse of PlantModules: `PlantSystem`.
"""

# ╔═╡ a3c5dba8-8acc-424a-87bc-d1c6a263578c
plantsys = PlantModules.PlantSystem(
	model_defaults = model_defaults,
	module_defaults = module_defaults,
	module_coupling = module_coupling,
	struct_connections = struct_connections,
	func_connections = func_connections
)

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
prob = ODEProblem(plantsys, time_span)

# ╔═╡ c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
sol = solve(prob)

# ╔═╡ a6608eff-9399-443c-a33a-c62341f7b14c
md"""
## Answering the toy problem

Finding the answer to our toy problem now comes down to plotting out the soil water content and visually inspecting when it gets too low:
"""

# ╔═╡ fe4df2d4-878e-41aa-8860-991c891e2dd2
plot(sol, struct_modules = [:Soil], func_vars = [:W]) #! imagine that this works

# ╔═╡ 2d131155-f62b-4f4a-92d8-9e7443202434
md"""
(Short discussion of plot and how it saved our lives)
"""

# ╔═╡ fb3c58df-1d6b-4ced-803d-2d0fc537b942
md"""
# TO-DO list

- Update hyperlinks with pages of docs when they exist
- Change parameter - and initial values with logical ones based on some data or something
- Add shapes in some better way
- Add Unitful units?
- How to add multiple soil compartments?
- Add connection between soil and air?
"""

# ╔═╡ 38c69eea-a4dd-4fc0-951f-dc36e9530b80


# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─d16a6d55-1f29-4b98-b1af-2dee1d38f386
# ╟─6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
# ╠═1144887a-a4c7-46f6-9cf8-cba50cf873d0
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
# ╟─4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
# ╟─3035b6d0-bca0-4803-b32a-da1459bdd880
# ╟─df4cd3de-d2b2-4f11-b755-a36e640fd2d5
# ╠═1394e5ed-39f3-4a9f-8737-be7183219314
# ╟─fbe00b62-6ca7-4c90-9050-081312911c74
# ╠═271d48a7-7022-4766-83d9-a70fab92515e
# ╟─3c13c600-5135-44ea-8fc2-a1e11f72d0c5
# ╠═5f21a4b0-f663-4777-94f3-d00acba109b3
# ╟─d574346b-e8d8-49cd-a3bc-bf53617119f1
# ╟─930e7ed8-0bfe-4e5a-8890-a1d1ce155881
# ╟─4cedbd9d-84ed-46f3-9a10-6cb993643f87
# ╠═d54705b3-d8f4-4cc2-a780-369343749113
# ╟─5885b970-cfb6-4e3b-8039-23a5f7ab4bae
# ╟─b0b5e539-19f5-4feb-8a2f-a6b2f0590c7c
# ╟─c4d42004-20dd-40ec-9200-1dd2b69b0bec
# ╠═9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
# ╟─ecb23f1e-ee39-4c5a-911f-eaa252c23968
# ╠═86d81fa3-bea4-40fa-9522-7db9fe2f6a82
# ╟─c36668fa-3046-4967-a616-841a048ea7f9
# ╠═3bf42137-1551-44d6-b7af-eab13a97b6ef
# ╠═db8fe96d-c8c2-47c0-9377-281ce577a920
# ╠═f267ca25-5b7f-48e4-9865-1c7952215a78
# ╠═61bf737a-2226-42dc-b93a-a8f578048268
# ╟─668e02ee-ac78-4b3d-983d-402ec04584ef
# ╠═caab574a-05c5-4c0d-9ae4-19fd514a6c6c
# ╟─f03a61ce-a0ff-43ff-abdd-2342f76e8c93
# ╟─ac6e098d-98ff-4940-8386-dc76c75eb2c3
# ╠═611289e9-e22c-4e6e-beec-ccea90eb48c9
# ╟─fb72735b-3d45-438d-ad83-3e36f42f5bb8
# ╟─210d81ef-153e-4744-8266-80af4099770c
# ╟─bc7573e7-bcd6-4347-9b0c-9111a824c9b5
# ╠═a3c5dba8-8acc-424a-87bc-d1c6a263578c
# ╟─d51795b2-32d3-455c-b785-5e501cfbdd08
# ╟─d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
# ╟─6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
# ╠═bf114636-1e35-49f1-9407-f472b443a9ea
# ╠═50d6fc31-80f5-4db7-b716-b26765008a0d
# ╠═c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
# ╟─a6608eff-9399-443c-a33a-c62341f7b14c
# ╠═fe4df2d4-878e-41aa-8860-991c891e2dd2
# ╟─2d131155-f62b-4f4a-92d8-9e7443202434
# ╟─fb3c58df-1d6b-4ced-803d-2d0fc537b942
# ╠═38c69eea-a4dd-4fc0-951f-dc36e9530b80
