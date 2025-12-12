# # Tutorial 1: Package basics

# ## Introduction

# We'll introduce the package by modeling the growth of a small pepper seedling. The goal is to simulate water transport to estimate when the plant will need watering again.
using Plots
using Revise, Infiltrator
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, OrdinaryDiffEq

# pretty colors :)
include(homedir() * raw"\Documents\Github\Caverns_of_code\Julia\Lifehacks\catpuccin\get_palette.jl")
cpalette = get_palette("prettycolors")

# ## Structure

# ### The plant

# The plant structure is defined as a **graph**. For this tutorial, we'll make use of [PlantGraphs.jl](https://github.com/VirtualPlantLab/PlantGraphs.jl)'s graph implementation.
# The different types of nodes in the graph represent different (repeating) structural parts of the plant, which we here refer to as **structural modules**.

struct Stem <: Node end
struct Leaf <: Node
	D::Vector
end

# Individual graph nodes can contain parameter values and initial values that should be used for their functional modules.
# Here, we'll give the leaves a size field `D` so we can start off one of them larger than the other.

plant_graph = Stem() + Stem() +
	(Leaf([5.0, 3.0, 0.05]), Leaf([1.0, 0.5, 0.05]));

plotstructure(plant_graph)

# ### The environment

# Likewise, the environment can also be modelled as graphs. This can also be included in the graph containing the plant structure.
# However, not all graph implementations (such as PlantGraphs.jl) support cyclical graphs. 
# This is why PlantModules allows entering multiple (acyclic) graphs and specifying how those graphs are connected, thereby allowing for cycles in the final structure.

struct Soil <: Node end
struct Air <: Node end

graphs = [Soil(), plant_graph, Air()]

# The connections between graphs are defined as a pair between the indices of the two graphs in question and what nodes to connect.
# The latter can be either the node itself, the structural module of the node, or a function that takes two nodes and returns whether they should be connected.
intergraph_connections = [[1, 2] => (:Soil, getnodes(plant_graph)[1]), [2, 3] => (:Leaf, :Air)];

# Finally we can combine all structural information in one variable:
plantstructure = PlantStructure(graphs, intergraph_connections)

begin
	using Random
	Random.seed!(12)
	plotstructure(plantstructure, palette = cpalette)
	# savefig(homedir() * raw"\Documents\Github\Den_of_evil\Non-note files\images" * "/fig_plantmodules_ex1_structure.pdf")
end

# ## Function

# ### Defining functionality

# Now that plant's structural modules are defined, we need to define their functional behaviour.
# This is done in the form of [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) `System`s.
# You can define new ones yourself, but here we'll use ones already defined by PlantModules. You can consult a list of available functional modules [here](https://www.youtube.com/watch?v=xvFZjo5PgG0).

# ### Coupling
# Our model needs to know which structural modules make use of which functional modules:
module_coupling = Dict(
	:Stem => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, simple_photosynthesis_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module]
);

# Aside from the graph's nodes, the edges also need to be assigned functionality.
# For this model, that comes down to describing how water flows from one plant part to another.
# We can use one of PlantModules' predefined functional modules again as follows:
connecting_modules = Dict(
	(:Soil, :Stem) => hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Leaf) => hydraulic_connection,
	(:Leaf, :Air) => daynight_hydraulic_connection,
);

plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# ### Parameter and initial values

# The parameters and initial values of the ODE systems defining the different plant parts are generally not identical.
# We can therefore change them as shown below.

# When multiple specifications for a value exist, PlantModules follows the following list of priority:
# > **Model-wide default values / module-specific default values / node-specific values**

# For __every node__:
default_changes = Dict(:Ψ => -0.2);

# For every node of a given __structural module__:
module_defaults = Dict(
	:Stem => Dict(:D => [0.5, 10.0], :M => 200e-6),
	:Leaf => Dict(:shape => Cuboid()),
	:Soil => Dict(:W_max => 1e3, :K => 1.0),
	:Air => Dict(:W_r => 0.7, :K => 1e-3)
);

# Note that changing the values for specific nodes is also possible, as discussed during the section on graph creation. 

# The parameters of the connections can also be changed.
connection_values = Dict(
	(:Leaf, :Air) => Dict(:t_sunrise => 6, :t_sunset => 22),
);

# All functional information also gets bundled, though the constructor is more complex.
# Check out its help page to see all the possible arguments and their defaults. 
plantparams = PlantParameters(; module_defaults, default_changes, connection_values)

# ## Generate and run system

# All that's left to do is running `generate_system` to get the plant's ODESystem
system = generate_system(plantstructure, plantcoupling, plantparams)

# ...and solving it.
time_span = (0, 48.0);
prob = ODEProblem(system, [], time_span, sparse = true, use_scc = false)
@time sol = solve(prob);

air_W_r = get_subsystem_variables(system, plantstructure, :W_r, :Air)[1]
prob2 = remake(prob, u0 = [air_W_r => 0.5])
sol2 = solve(prob2)

# ## Plotting
# Finally, we can use PlantModules' `plotgraph` function to more easily plot the desired results.
# Based on the water content `W` of the soil (which was the second graph), and growth of the leaves, we can plan when to water next!
p1 = Plots.plot(
	plotgraph(sol, plantstructure, varname = :ΣF, structmod = :Air, title = "Low transpiration", xaxis = "t (hr)"),
	plotgraph(sol2, plantstructure, varname = :ΣF, structmod = :Air, title = "High transpiration", xaxis = "t (hr)"), ylims = (0.0, 0.25), yaxis = "ΣF"
)

p2 = Plots.plot(
	plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Stem, :Leaf], title = "Low transpiration"),
	plotgraph(sol2, plantstructure, varname = :Ψ, structmod = [:Stem, :Leaf], title = "High transpiration"),
	ylims = (-0.35, -0.1), yaxis = "Ψ"
)

begin
	Plots.plot(
		plotgraph(sol, plantstructure, varname = :ΣF, structmod = :Air, title = "Low evaporative demand", ylims = (0.0, 0.25), yaxis = "Net water influx (g / hr)", palette = cpalette, lw = 2),
		plotgraph(sol2, plantstructure, varname = :ΣF, structmod = :Air, title = "High evaporative demand", ylims = (0.0, 0.25), palette = cpalette, lw = 2),
		plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Stem, :Leaf], title = "", ylims = (-0.4, -0.1), yaxis = "Water potential (MPa)", xaxis = "Time (hr)", palette = cpalette, lw = 2),
		plotgraph(sol2, plantstructure, varname = :Ψ, structmod = [:Stem, :Leaf], title = "", ylims = (-0.4, -0.1), xaxis = "Time (hr)", palette = cpalette, lw = 2),
		size = (800, 600)
	)
	# savefig(homedir() * "\\Documents\\Github\\Den_of_evil\\Non-note files\\images\\" * "fig_plantmodules_ex1_results.pdf")
end