module PlantModules

# # Imports
using Graphs # used to define plant structure
import Graphs: edges, edgetype, has_edge, has_vertex, ne, nv, neighbors,
    inneighbors, outneighbors, vertices, is_directed, src, dst # To define own Graphs.jl `AbstractGraph`
using SciMLBase # interface with ODE solver packages
using ModelingToolkit, Unitful # simulation
import ModelingToolkit: get_eqs, get_unknowns, get_ps, get_parameter_dependencies, get_observed,
    get_continuous_events, get_discrete_events, get_defaults, get_systems, get_name, get_iv,
    get_gui_metadata, get_parent, get_systems # MTK internals
import ModelingToolkit: parameter_values, setp # MTK internals for problem remaking
import ModelingToolkit.SciMLBase.AbstractSciMLProblem
using PlantGraphs, MultiScaleTreeGraph # graphs
using RecipesBase, GraphRecipes # visualisation

# # Exports
# ## Own functions
export generate_system, PlantStructure, PlantCoupling, PlantParameters # generate system
export ModuleShape, Sphere, Cylinder, Cuboid # shapes
export getnodes, getneighbors, getattributes, getstructmod, getid # graph functions
export getdimensionality, volume, cross_area, surface_area # shape functions
export logsumexp, smooth_daynight # smooth functions
export hydraulic_module, environmental_module, constant_carbon_module, simple_photosynthesis_module, Ψ_soil_module, Ψ_air_module, K_module, constant_K_module # node modules
export hydraulic_connection, constant_hydraulic_connection, daynight_hydraulic_connection # edge modules
export readXEG, convert_to_MTG, convert_to_PG # graph reading and converting
export remake_graphsystem, remake_graphsystem!, get_subsystem_variables # system remaking
export plotstructure, plotgraph, plotnode # plotting

# ## Re-exports
export graphplot

# # Include src files
include("plantstructure.jl")
include("graph_functions.jl")
include("graph_reading.jl")
include("graph_conversion.jl")
include("generate_system.jl")
include("shapes.jl")
include("smoothfuncs.jl")
include("func_modules.jl")
include("system_remaking.jl")
include("plotting.jl")

end # module