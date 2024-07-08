### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 57b8dcb8-9baa-4ddf-9368-431b1be5850e
using Pkg; Pkg.activate(".")

# ╔═╡ 662c0476-70aa-4a60-a81c-d7db2248b728
include("../src/PlantModules.jl"); using .PlantModules

# ╔═╡ 65f88593-1180-447a-900f-49aef4647cd1
using PlantGraphs, ModelingToolkit, DifferentialEquations, Plots

# ╔═╡ 56c3527f-d8df-4f5c-9075-77c34d5c7204
md"""
# Tutorial 1: Package basics
"""

# ╔═╡ 6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
md"""
## Introduction 👋

In this tutorial, we'll cover the basics of using the PlantModules package. Here's a brief overview:
- Toy problem: when to water your pepper seedlings
- Creating a modular plant model
  - Defining structural and functional plant modules
  - Specifying model parameters
  - Describing how the modules interact
- Running the model
- Answering the toy question
"""

# ╔═╡ 1144887a-a4c7-46f6-9cf8-cba50cf873d0
md"""
### Toy problem description

There is no fun in modeling without purpose, which is why we'll introduce the package with a simple problem.

Our recently sprouted pepper seedlings are growing on the windowsill indoors. We just watered them, and we'd like to know how the water content of the soil changes through time so we can water them again at the optimal moment: the soil should have dried adequately to promote root growth without putting the seedlings under too much drought stress.

![plantfigu](https://www.almanac.com/sites/default/files/users/The%20Editors/pepper-seedlings_full_width.jpg)
"""
#! replace picture with something royalty free. a drawing?

# ╔═╡ bc44fae7-6c4a-4f27-be51-dbabf4037601
md"""
### Loading the necessary packages

Before we can get started, we need to do some basic setup: activating our project and loading in the required packages.
"""

# ╔═╡ c2970976-cd1d-4e89-9df2-83932f934f2d
import GLMakie.draw

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

# ╔═╡ 0cc02e82-4fe8-4f27-a2d2-eb4bfba6b291
md"""
Then we can define the structural modules of our plant. For our example, we will define three structural modules on the organ scale of the plant.
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
There are five functional modules already provided by PlantModules:

- `PlantModules.hydraulic_module`, which describes the hydraulics-driven growth of a compartment.
- `PlantModules.environmental_module`, which describes a compartment that is subject to hydraulic laws but does therefor not change in size.
- `PlantModules.constant_carbon_module`, which describes a compartment with a constant amount of carbon.
- `PlantModules.Ψ_soil_module`, which describes an empirical relationship between the total water potential and relative water content of soil.
- `PlantModules.Ψ_air_module`, which describes the relationship between the total water potential and water potential of air.
"""

# ╔═╡ 1e2eaa86-d6e0-4749-bc92-04c64fe9f47d
md"""
The first module contains the most important (and theoretically founded) functionality of the package, whereas the other four are mostly to allow for implementing a model that works out of the box. Users that want to create a realistic model are expected to implement some functional modules of their own in accordance with their own use for the package, as will be discussed next tutorial. For now, we will stick to the pre-implemented functional modules.
"""

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
PlantModules.default_params |> print # Parameter values

# ╔═╡ ae8b0cb6-6f0f-4c18-b05f-01aec542037c
PlantModules.default_u0s |> print # Initial values

# ╔═╡ fbe00b62-6ca7-4c90-9050-081312911c74
md"""
Now imagine we have some data on our plant in question and it calls for different default values. We could manually rewrite the defaults from scratch, but using the `alter_defaults` function is a lot faster: 
"""

# ╔═╡ 271d48a7-7022-4766-83d9-a70fab92515e
default_changes = (Γ = 0.4, P = 0.2, T = 293.15)

# ╔═╡ 2c516ca9-e169-43aa-b41b-b5f91ee80588
default_params, default_u0s = PlantModules.alter_defaults(default_changes)

# ╔═╡ e60a5281-18ce-4e88-b76b-b8279af7dc65
md"""
This changes the parameter- and initial values to the specified new value in every functional module:
"""

# ╔═╡ 8f5b2fb2-7409-4a25-b74f-bcfb247ae647
default_params |> print

# ╔═╡ 3c13c600-5135-44ea-8fc2-a1e11f72d0c5
md"""
Next to changing functional values over the entire model, some of the information we have may require different default values specifically for certain structural modules. We can specify this as shown below.
"""

# ╔═╡ de3220a3-18d9-423c-829f-2b43cac301fe
begin
	C_root = 300e-6
	C_stem = 400e-6
	C_leaf = 450e-6 
end

# ╔═╡ 5f21a4b0-f663-4777-94f3-d00acba109b3
module_defaults = (
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 10.0, 20.0], ϕ_D = [0.7, 0.1, 0.05]), D = [30, 3, 1], M = C_root),
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 15.0], ϕ_D = [0.8, 0.03]), D = [1.5, 10], M = C_stem),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 500.0, T = 288.15),
	Air = (W_r = 0.8,)
)

# ╔═╡ e0e5e7f7-d2f5-404b-a6c6-6848c318eccb
md"""
As you can see, we changed the default $shape$ between the root, stem and leaf modules, as well as their initial dimensions $D$ and metabolite concentration $M$. For our soil module, we changed the maximum water content and temperature. For a detailed explanation of the parameters and initial values here, we again refer to the [theoretical overview](nothinghere).
"""

# ╔═╡ 930e7ed8-0bfe-4e5a-8890-a1d1ce155881
md"""
### Coupling functional and structural modules
"""

# ╔═╡ 4cedbd9d-84ed-46f3-9a10-6cb993643f87
md"""
Our model needs to know which structural modules make use of which functional modules. As perhaps the simplest part of our modeling workflow, we can define this using a `Vector` of `Pairs`.
"""

# ╔═╡ d54705b3-d8f4-4cc2-a780-369343749113
module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
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

Since we're using the PlantGraphs Julia package for this tutorial, we can easily define simple graphs using its syntax:
"""

# ╔═╡ 9af27c17-8f21-4f22-a5bb-e9c95cfdf2f9
plant_graph = Root() + Stem() + (Leaf([2]), Leaf([2.5])) # Leaves are chosen to be  modelled as spheres and are here given radii of 2 and 2.5 cms

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

# ╔═╡ 8980f1eb-e461-4624-9cce-83a7cf822349
graphs = [plant_graph, soil_graph, air_graph]

# ╔═╡ 20049311-d6e6-41d3-a0d8-8bad88c174f9
md"""
The following syntax means: "connect all Root nodes from `graphs[1]` to all Soil nodes from `graphs[2]`", and so on. Note that the order matters! For example, `[1, 2] => (:Root, :Soil)` is not the same as `[2, 1] => (:Root, :Soil)`.
"""

# ╔═╡ 61bf737a-2226-42dc-b93a-a8f578048268
intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] # Let's also add a connection between the soil and the air to simulate direct evaporation

# ╔═╡ 668e02ee-ac78-4b3d-983d-402ec04584ef
md"""
The full structural connection information can then be acquired by putting them together in a vector.
"""

# ╔═╡ caab574a-05c5-4c0d-9ae4-19fd514a6c6c
struct_connections = [graphs, intergraph_connections]

# ╔═╡ f03a61ce-a0ff-43ff-abdd-2342f76e8c93
md"""
> For our tiny example plant, we only have one structural module that actually repeats, somewhat defeating the purpose of modelling in a modular manner: we may as well write this entire model out by hand! Rest assured, however, that the approach we're seeing here also works for larger plants, as we will see in the next tutorial.
"""

# ╔═╡ ac6e098d-98ff-4940-8386-dc76c75eb2c3
md"""
#### Functional connections

If we consider functional modules as the information stored in the nodes of the graph, then the functional connections are the information stored in the edges. For the functional process under consideration here, this corresponds with a ModelingToolkit equation set describing water flow between compartments. It only has one parameter *K*, which is the hydraulic conductivity of the connection in question. 

Defining a new functional connection module is discussed in the following tutorial, for now we'll simply change the parameter value of this one:
"""

# ╔═╡ 611289e9-e22c-4e6e-beec-ccea90eb48c9
connecting_modules = [
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 10]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Stem, :Leaf) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-2]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-2])
]

# ╔═╡ a8a725d4-d876-4867-acbc-26bbadc4b462
md"""
These connecting modules define the functional information of the graph edges. In this case: the flow of water between nodes based on some hydraulic conductivity parameter and a difference in water potentials.
However, the system also needs to know how the edges relate to the nodes functionally. For example, the net amount of water coming into a node is the sum of all water flows of edges connected to this node. This information is given under the form of a function returning a set of equations. For the hydraulic model, this is already provided.
"""

# ╔═╡ 5ca7ded4-d333-4edb-96db-3fdb7bc827ce
get_connection_eqs = PlantModules.multi_connection_eqs

# ╔═╡ 113c0bdc-53e4-4a19-a47c-f4afba237eeb
md"""
The functional connections are then defined as the combination of the functional modules and the function generating connecting equations.
"""

# ╔═╡ ebf46963-bb9d-4604-924c-2b051189debc
func_connections = [connecting_modules, get_connection_eqs]

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
Now that we have all parts of our model defined, all that's left to do is putting it together. For this we use the main workhorse of PlantModules: `generate_system`.
"""

# ╔═╡ a3c5dba8-8acc-424a-87bc-d1c6a263578c
system = PlantModules.generate_system(default_params, default_u0s, module_defaults, module_coupling, struct_connections, func_connections);

# ╔═╡ d51795b2-32d3-455c-b785-5e501cfbdd08
md"""
This function will generate the `ODESystem` describing the model. It is possible to fine-tune the model even further at this stage as described in the [Customizing the model](nothinghere) section of the docs, thought this should generally not be required.
"""

# ╔═╡ d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
md"""
## Running the model 🏃‍♂️
"""

# ╔═╡ 6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
md"""
The rest of the modeling workflow is mostly taken care of by ModelingToolkit.jl and DifferentialEquations.jl, with some syntactic sugar added by PlantModules. For users that are unfamiliar with the package, it is recommended to take a brief look at [the ModelingToolkit docs](https://docs.sciml.ai/ModelingToolkit/stable/) before proceeding.
"""

# ╔═╡ bf114636-1e35-49f1-9407-f472b443a9ea
time_span = (0, 7*24.0) # We'll simulate our problem for a timespan of one week

# ╔═╡ 2f431e8c-d0e4-4117-896f-3140d9633d1d
sys_simpl = structural_simplify(system);

# ╔═╡ 50d6fc31-80f5-4db7-b716-b26765008a0d
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), time_span);

# ╔═╡ c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
sol = solve(prob)

# ╔═╡ a6608eff-9399-443c-a33a-c62341f7b14c
md"""
## Answering the toy problem

Finding the answer to our toy problem now comes down to plotting out the soil water content and visually inspecting when it gets too low:
"""

# ╔═╡ e890700e-80a4-4dfc-8380-732bf91d1aa4
PlantModules.plotnode(sol, graphs[2], func_varname = :W_r)

# ╔═╡ 2d131155-f62b-4f4a-92d8-9e7443202434
md"""
It looks like after 1 week only about 20% of the soil water content remains. Since we don't want our seedlings to experience too much drought stress, we may want to water it when the relative soil water content is at about 50%, which would be after around half a week.

Of course, the functional modules used in this tutorial are a gross oversimplification of reality and these results are very iffy at best. Next tutorial we'll see how we can create our own functional processes to make the model more realistic.
"""

# ╔═╡ Cell order:
# ╟─56c3527f-d8df-4f5c-9075-77c34d5c7204
# ╟─6ab177fd-ed5b-4ae4-a2b5-f7f4eb8e4d0d
# ╟─1144887a-a4c7-46f6-9cf8-cba50cf873d0
# ╟─bc44fae7-6c4a-4f27-be51-dbabf4037601
# ╠═57b8dcb8-9baa-4ddf-9368-431b1be5850e
# ╠═662c0476-70aa-4a60-a81c-d7db2248b728
# ╠═65f88593-1180-447a-900f-49aef4647cd1
# ╠═c2970976-cd1d-4e89-9df2-83932f934f2d
# ╟─aa3b75e4-1868-4c84-8dc8-f9b54d560b3a
# ╟─6ef5c63a-b753-43ae-baee-f6c24313a385
# ╟─b6eb66b5-a2d7-4baf-b6a6-87e819309a2d
# ╟─aec7bcd6-6f27-4cf5-a955-f4d59e778fd3
# ╟─659e911c-8af2-4a66-855a-e333c41120c1
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
# ╟─1e2eaa86-d6e0-4749-bc92-04c64fe9f47d
# ╟─4d17b269-06b8-4293-b2cb-b6bd9fa0ccc8
# ╟─3035b6d0-bca0-4803-b32a-da1459bdd880
# ╟─df4cd3de-d2b2-4f11-b755-a36e640fd2d5
# ╠═1394e5ed-39f3-4a9f-8737-be7183219314
# ╠═ae8b0cb6-6f0f-4c18-b05f-01aec542037c
# ╟─fbe00b62-6ca7-4c90-9050-081312911c74
# ╠═271d48a7-7022-4766-83d9-a70fab92515e
# ╠═2c516ca9-e169-43aa-b41b-b5f91ee80588
# ╟─e60a5281-18ce-4e88-b76b-b8279af7dc65
# ╠═8f5b2fb2-7409-4a25-b74f-bcfb247ae647
# ╟─3c13c600-5135-44ea-8fc2-a1e11f72d0c5
# ╠═de3220a3-18d9-423c-829f-2b43cac301fe
# ╠═5f21a4b0-f663-4777-94f3-d00acba109b3
# ╟─e0e5e7f7-d2f5-404b-a6c6-6848c318eccb
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
# ╠═8980f1eb-e461-4624-9cce-83a7cf822349
# ╟─20049311-d6e6-41d3-a0d8-8bad88c174f9
# ╠═61bf737a-2226-42dc-b93a-a8f578048268
# ╟─668e02ee-ac78-4b3d-983d-402ec04584ef
# ╠═caab574a-05c5-4c0d-9ae4-19fd514a6c6c
# ╟─f03a61ce-a0ff-43ff-abdd-2342f76e8c93
# ╟─ac6e098d-98ff-4940-8386-dc76c75eb2c3
# ╠═611289e9-e22c-4e6e-beec-ccea90eb48c9
# ╟─a8a725d4-d876-4867-acbc-26bbadc4b462
# ╠═5ca7ded4-d333-4edb-96db-3fdb7bc827ce
# ╟─113c0bdc-53e4-4a19-a47c-f4afba237eeb
# ╠═ebf46963-bb9d-4604-924c-2b051189debc
# ╟─fb72735b-3d45-438d-ad83-3e36f42f5bb8
# ╟─210d81ef-153e-4744-8266-80af4099770c
# ╟─bc7573e7-bcd6-4347-9b0c-9111a824c9b5
# ╠═a3c5dba8-8acc-424a-87bc-d1c6a263578c
# ╟─d51795b2-32d3-455c-b785-5e501cfbdd08
# ╟─d3d7b52b-016b-4c17-a4cc-18ec4ad8d686
# ╟─6b46bf1d-b54e-48e3-b4eb-364b4e2b1dfd
# ╠═bf114636-1e35-49f1-9407-f472b443a9ea
# ╠═2f431e8c-d0e4-4117-896f-3140d9633d1d
# ╠═50d6fc31-80f5-4db7-b716-b26765008a0d
# ╠═c38b1a71-c5e9-4bfa-a210-bcbf9068f7ed
# ╟─a6608eff-9399-443c-a33a-c62341f7b14c
# ╠═e890700e-80a4-4dfc-8380-732bf91d1aa4
# ╟─2d131155-f62b-4f4a-92d8-9e7443202434
