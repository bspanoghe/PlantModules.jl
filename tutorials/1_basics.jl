# # Tutorial 1: Package basics

# ## Introduction

# We'll introduce the package by modeling the growth of a small pepper seedling. The goal is to simulate water transport to estimate when the plant will need watering again.

using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Plots
import GLMakie.draw

# ## Structure

# ### The plant

# The plant structure is defined as a **graph**. For this tutorial, we'll make use of [PlantGraphs.jl](https://github.com/VirtualPlantLab/PlantGraphs.jl)'s graph implementation.
# The different types of nodes in the graph represent different (repeating) structural parts of the plant, which we here refer to as **structural modules**.

mutable struct Root <: Node end
mutable struct Stem <: Node end
mutable struct Leaf <: Node
	D::Vector
end

# Individual graph nodes can contain parameter values and initial values that should be used for their functional modules.
# Here, we'll give the leaves a size field `D` so we can start off one of them larger than the other.

plant_graph = Root() + Stem() + (Leaf([5, 2, 0.01]), Leaf([4, 1.5, 0.01]))
draw(plant_graph, resolution = (500, 400))

# ### The environment

# Likewise, the environment can also be modelled as graphs. This can also be included in the graph containing the plant structure.
# However, not all graph implementations (such as PlantGraphs.jl) support cyclical graphs. 
# This is why PlantModules allows entering multiple (acyclic) graphs and specifying how those graphs are connected, thereby allowing for cycles in the final structure.

struct Soil <: Node end
struct Air <: Node end

graphs = [plant_graph, Soil(), Air()]

# The connections between graphs are defined as a pair between the indices of the two graphs in question and what nodes to connect.
# The latter can be either the node itself, the structural module of the node, or a function that takes two nodes and returns whether they should be connected.
intergraph_connections = [[1, 2] => (PlantModules.root(plant_graph), :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)]

# Finally we can combine all structural information in one variable:
struct_connections = PlantStructure(graphs, intergraph_connections)


# ## Function

# ### Defining functionality

# Now that plant's structural modules are defined, we need to define their functional behaviour.
# This is done in the form of [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) `ODESystem`s.
# You can define new ones yourself, but here we'll use ones already defined by PlantModules. You can consult a list of available functional modules [here](https://www.youtube.com/watch?v=xvFZjo5PgG0).

# ### Finetuning

# The parameters and initial values of the ODE systems defining the different plant parts are generally not identical.
# We can therefor change them as shown below.

# For __every node__:
model_wide_changes = Dict(:Γ => 0.4, :P => 0.2, :T => 293.15) 

# For every node of a given __structural module__:
module_defaults = Dict(
	:Root => Dict(:M => 300e-6),
	:Stem => Dict(:D => [1.5, 10], :M => 400e-6),
	:Leaf => Dict(:shape => Cuboid(ϵ_D = [5.0, 10.0, 50.0], ϕ_D = [5e-3, 1e-3, 5e-5]), :M => 450e-6),
	:Soil => Dict(:W_max => 500.0, :T => 288.15, :W_r => 1.0),
) 

# Changing the values for specific nodes is also possible, as discussed during the section on graph creation. 

# When multiple specifications for a value exist, PlantModules follows the following list of priority:
# > **Model-wide default values / module-specific default values / node-specific values**

# ### Edge functionality

# Aside from the graph's nodes, the edges also need to be assigned functionality.
# For this model, that comes down to describing how water flows from one plant part to another.
# We can use one of PlantModules' predefined functional modules again and simply change the parameters as follows:
connecting_modules = [
	(:Soil, :Root) => (const_hydraulic_connection, Dict(:K => 10)),
	(:Root, :Stem) => (const_hydraulic_connection, Dict(:K => 5)),
	(:Stem, :Leaf) => (const_hydraulic_connection, Dict(:K => 5)),
	(:Leaf, :Air) => (const_hydraulic_connection,  Dict(:K => 2e-2)),
	(:Soil, :Air) => (const_hydraulic_connection,  Dict(:K => 5e-3))
]

# All functional information also gets bundled, though the constructor is more complex.
# Check out its help page to see all the possible arguments and their defaults. 
func_connections = PlantFunctionality(module_defaults = module_defaults,
	connecting_modules = connecting_modules, default_changes = model_wide_changes
)

# ## Coupling
# As the final piece of required information, our model needs to know which structural modules make use of which functional modules:
module_coupling = Dict(
	:Root => [hydraulic_module, constant_carbon_module],
	:Stem => [hydraulic_module, constant_carbon_module],
	:Leaf => [hydraulic_module, constant_carbon_module],
	:Soil => [environmental_module, Ψ_soil_module],
	:Air => [environmental_module, Ψ_air_module],
)

# ## Generate and run system

# All that's left to do is running `generate_system` to get the plant's ODESystem
system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false)

# ...and solving it.
# We refer to [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)'s documentation for more information on this part.
sys_simpl = structural_simplify(system)
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
sol = solve(prob);

# ## Plotting
# Finally, we can use PlantModules' `plotgraph` function to more easily plot the desired results.
# Based on the water content `W` of the soil (which was the second graph), we can plan when to water next!
plotgraph(sol, graphs[2], varname = :W)

# Some other variables we may be interested in, to showcase the plotting functionality:
# the growth of the leaves
plotgraph(sol, graphs[1], varname = :V, structmod = :Leaf)
# the water potentials of the plant parts and the soil
plotgraph(sol, graphs[1:2], varname = :Ψ)
# and the net water flux in all components of the system
plotgraph(sol, graphs, varname = :ΣF)