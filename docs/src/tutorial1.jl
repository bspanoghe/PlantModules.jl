# # Tutorial 1: Package basics

# ## Introduction 👋

# In this tutorial, we'll cover the basics of using the PlantModules package. Here's a brief overview:
# - Toy problem: when to water your pepper seedlings
# - Creating a modular plant model
#   - Defining structural and functional plant modules
#   - Specifying model parameters
#   - Describing how the modules interact
# - Running the model
# - Answering the toy question


# ### Toy problem description

# There is no fun in modeling without purpose, which is why we'll introduce the package with a simple problem.

# Our recently sprouted pepper seedlings are growing on the windowsill indoors. We just watered them, and we'd like to know how the water content of the soil changes through time so we can water them again at the optimal moment: the soil should have dried adequately to promote root growth without putting the seedlings under too much drought stress.

# ### Loading the necessary packages

# Before we can get started, we need to do some basic setup: activating our project and loading in the required packages.

using Pkg; Pkg.activate(".")
using PlantModules
using PlantGraphs, ModelingToolkit, DifferentialEquations, Plots

import GLMakie.draw

# ## Creating the model 🛠


# As the name suggests, one of the most important goals of PlantModules is to enable the user to easily model plant growth in a _modular_ manner. In order to achieve this, we will define our plant in function of a few sets of similarly behaving parts or "modules" for short. Parts of the same module can be similar on the structural - and on the functional level.

# On the structural level, some examples are:
# - An oak tree can be considered a repeating module in a forest.
# - A branch can be considered a repeating module in a tree.
# - A collenchyma cell can be considered a repeating module in a branch.

# On the functional level, on the other hand:
# - Most structural modules of a tree (branches, leaves, root segments, etc.) contain water which will flow between them as dictated by hydraulic laws, which plays an important role in plant growth. All these structural modules share this same functional module. 
# - One or more of the plant's structural modules will assimilate carbon through photosynthesis, while the others do not. On the cellular level, only cells containing chloroplasts should share this functional module.

# As such, the workflow to create a model in PlantModules boils down to defining these modules and how they interact. Let's jump right in!

# ### Defining the structural modules
# #### The plant

# The first step in modeling our plant's growth is defining its structure. There are a lot of options for formalizing the structure of a plant, and perhaps the most obvious first choice we need to make is what spatial scale to model on. PlantModules allows for a lot of freedom here, including the option to combine multiple spatial scales. 
# For the first tutorial, we'll consider a very simple example: an organ-scale plant model of a plant with two leaves, a stem and a root system.
# PlantModules uses **graphs** to define the relationships between modules. As such, structural modules need to be implemented as graph nodes of some sort. We'll be using the graph implementation from [the PlantGraphs](https://virtualplantlab.com/stable/manual/Graphs/) package for this tutorial, though any graph implementation can be used.
# Then we can define the structural modules of our plant. For our example, we will define three structural modules on the organ scale of the plant.

mutable struct Root <: Node end
mutable struct Stem <: Node end
mutable struct Leaf <: Node
	D::Vector
end

# Individual graph nodes can contain parameter values and initial values to allow differences between nodes of the same type. Here, we'll give the leaves a size field `D` so we can start off one of them larger than the other. We'll see more on parameter - and initial values later.
# #### The environment

# Plants need an environment to grow in. For most plants, the most basic environmental compartments that need defining are the soil, from which the plant gets water and nutrients, and the air, with which the plant exchanges gasses and to which it loses water. 

# For more intricate models one may for example also want to take into account the sun as well as possible sources of shadow or divide the soil into multiple compartments, but we will stick to the basics for now.

struct Soil <: Node end
struct Air <: Node end

# ### Defining the functional modules

# Now that we have an idea of the plant's structural modules, we need to define some sort of functionality expressed by one or more of them. Some examples, such as only a select set of structural modules doing photosynthesis, have been mentioned before. In this tutorial, we'll mostly neglect carbon dynamics and instead focus on the plant's water dynamics.
# PlantModules defines functional modules as as sets of differential equations implemented in [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/).
# There are five functional modules already provided by PlantModules:

# - `PlantModules.hydraulic_module`, which describes the hydraulics-driven growth of a compartment.
# - `PlantModules.environmental_module`, which describes a compartment that is subject to hydraulic laws but does therefor not change in size.
# - `PlantModules.constant_carbon_module`, which describes a compartment with a constant amount of carbon.
# - `PlantModules.Ψ_soil_module`, which describes an empirical relationship between the total water potential and relative water content of soil.
# - `PlantModules.Ψ_air_module`, which describes the relationship between the total water potential and water potential of air.
# The first module contains the most important (and theoretically founded) functionality of the package, whereas the other four are mostly to allow for implementing a model that works out of the box. Users that want to create a realistic model are expected to implement some functional modules of their own in accordance with their own use for the package, as will be discussed next tutorial. For now, we will stick to the pre-implemented functional modules.

# #### Defining parameter and initial values 
# Next to the equations themselves, an important part in describing the plant's functional processes is defining the equations' parameter - and initial values. Since we'll restrict ourselves to hydraulic equations for this tutorial, these parameters and variables are the ones described in the section [Theoretical overview](@ref).

# In PlantModules, parameter - and initial values are defined in a hierarchical manner, here listed in ascending priority:

# > **Model-wide default values < module-specific default values < node-specific values**

# For our current example, this means that there is no point in specifying the initial values for the compartment dimensions *D* for our Leaf compartments in the module-specific default values or the model-wide default values, since we already defined these values in all of the Leaf nodes of the graph, which have the highest priority.
# Before we change any of the default parameter - and initial values, we can take a look at the defaults used by PlantModules:

PlantModules.default_params |> print # Parameter values
PlantModules.default_u0s |> print # Initial values

# Now imagine we have some data on our plant in question and it calls for different default values. We could manually rewrite the defaults from scratch, but using the `alter_defaults` function is a lot faster: 

default_changes = (Γ = 0.4, P = 0.2, T = 293.15)
default_params, default_u0s = PlantModules.alter_defaults(default_changes)

# This changes the parameter- and initial values to the specified new value in every functional module:

default_params |> print

# Next to changing functional values over the entire model, some of the information we have may require different default values specifically for certain structural modules. We can specify this as shown below.

C_root = 300e-6
C_stem = 400e-6
C_leaf = 450e-6 

module_defaults = (
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 10.0, 20.0], ϕ_D = [0.7, 0.1, 0.05]), D = [30, 3, 1], M = C_root),
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 15.0], ϕ_D = [0.8, 0.03]), D = [1.5, 10], M = C_stem),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 500.0, T = 288.15),
	Air = (W_r = 0.8,)
)

# As you can see, we changed the default $shape$ between the root, stem and leaf modules, as well as their initial dimensions $D$ and metabolite concentration $M$. For our soil module, we changed the maximum water content and temperature. For a detailed explanation of the parameters and initial values here, we again refer to the [Theoretical overview](theory.md).

# ### Coupling functional and structural modules

# Our model needs to know which structural modules make use of which functional modules. As perhaps the simplest part of our modeling workflow, we can define this using a `Vector` of `Pairs`.

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

# ### Module connections
# The final information required to fully define our plant growth model is how the modules relate to each other. This can once again be separated into a structural and functional part.

# #### Structural connections

# As mentioned near the beginning, PlantModules uses graphs to define relationships between modules, and the structural modules were implemented as graph node types. All we need to do now is create a graph using these node types reflecting the desired structure of our plant.

# Since we're using the PlantGraphs Julia package for this tutorial, we can easily define simple graphs using its syntax:

plant_graph = Root() + Stem() + (Leaf([2]), Leaf([2.5])) # Leaves are chosen to be  modelled as spheres and are here given radii of 2 and 2.5 cms

# We can inspect the structure of our plant by visualising the graph:

draw(plant_graph, resolution = (500, 400))

# There are some structural modules missing here: the environmental ones. We could have added the soil compartment in the graph we just defined, but not the air. That is because both leaves are connected to this same air node, which would introduce a cycle into our graph, something many plant-related graph implementations do not support.

# To overcome this, PlantModules allows the user to input multiple graphs and the ability to define what structural modules should be fully connected between all input graphs to acquire the final structural graph:

soil_graph = Soil()
air_graph = Air()
graphs = [plant_graph, soil_graph, air_graph]

# The following syntax means: "connect all Root nodes from `graphs[1]` to all Soil nodes from `graphs[2]`", and so on. Note that the order matters! For example, `[1, 2] => (:Root, :Soil)` is not the same as `[2, 1] => (:Root, :Soil)`.

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] # Let's also add a connection between the soil and the air to simulate direct evaporation

# The full structural connection information can then be acquired by putting them together in a vector.

struct_connections = [graphs, intergraph_connections]
# > For our tiny example plant, we only have one structural module that actually repeats, somewhat defeating the purpose of modelling in a modular manner: we may as well write this entire model out by hand! Rest assured, however, that the approach we're seeing here also works for larger plants, as we will see in the next tutorial.

# #### Functional connections

# If we consider functional modules as the information stored in the nodes of the graph, then the functional connections are the information stored in the edges. For the functional process under consideration here, this corresponds with a ModelingToolkit equation set describing water flow between compartments. It only has one parameter *K*, which is the hydraulic conductivity of the connection in question. 

# Defining a new functional connection module is discussed in the following tutorial, for now we'll simply change the parameter value of this one:

connecting_modules = [
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 10]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Stem, :Leaf) => (PlantModules.hydraulic_connection, [:K => 5]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-2]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-2])
]

# These connecting modules define the functional information of the graph edges. In this case: the flow of water between nodes based on some hydraulic conductivity parameter and a difference in water potentials.
# However, the system also needs to know how the edges relate to the nodes functionally. For example, the net amount of water coming into a node is the sum of all water flows of edges connected to this node. This information is given under the form of a function returning a set of equations. For the hydraulic model, this is already provided.

get_connection_eqs = PlantModules.multi_connection_eqs

# The functional connections are then defined as the combination of the functional modules and the function generating connecting equations.

func_connections = [connecting_modules, get_connection_eqs]
# > ⚠ In contrast to the nodes, PlantModules currently does not support edges of the same type (that is, between the same structural modules) to specify different parameter - or initial values in the graph definition.

# ### Bringing it all together

# Now that we have all parts of our model defined, all that's left to do is putting it together. For this we use the main workhorse of PlantModules: `generate_system`.

system = PlantModules.generate_system(default_params, default_u0s, module_defaults, module_coupling, struct_connections, func_connections);

# This function will generate the `ODESystem` describing the model.

# ## Running the model 🏃‍♂️
# The rest of the modeling workflow is mostly taken care of by ModelingToolkit.jl and DifferentialEquations.jl, with some syntactic sugar added by PlantModules. For users that are unfamiliar with the package, it is recommended to take a brief look at [the ModelingToolkit docs](https://docs.sciml.ai/ModelingToolkit/stable/) before proceeding.
time_span = (0, 7*24.0) # We'll simulate our problem for a timespan of one week
sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), time_span);
sol = solve(prob)
# ## Answering the toy problem

# Finding the answer to our toy problem now comes down to plotting out the soil water content and visually inspecting when it gets too low:
PlantModules.plotnode(sol, graphs[2], func_varname = :W_r)
# It looks like after 1 week only about 20% of the soil water content remains. Since we don't want our seedlings to experience too much drought stress, we may want to water it when the relative soil water content is at about 50%, which would be after around half a week.

# Of course, the functional modules used in this tutorial are a gross oversimplification of reality and these results are very iffy at best. Next tutorial we'll see how we can create our own functional processes to make the model more realistic.