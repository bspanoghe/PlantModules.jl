module PlantModules

# # Exports
export generate_system, PlantStructure, PlantFunctionality # generate system
export Shape, Sphere, Cylinder, Cuboid # shapes
export getnodes, getroot # graph functions
export getdimensionality, volume, cross_area, surface_area # shape functions
export logsumexp, smooth_daynight # smooth functions
export hydraulic_module, environmental_module, constant_carbon_module, daynight_carbon_module, Ψ_soil_module, Ψ_air_module, K_module, constant_K_module # node modules
export hydraulic_connection, constant_hydraulic_connection, evaporation_connection # edge modules
export readXEG, convert_to_MTG, convert_to_PG # graph reading and converting
export plotgraph, plotnode # plotting

# # Imports
using SciMLBase # interface with ODE solver packages
using ModelingToolkit, Unitful # simulation
import ModelingToolkit: get_eqs, get_unknowns, get_ps, get_parameter_dependencies, get_observed,
    get_continuous_events, get_discrete_events, get_defaults, get_systems, get_name, get_iv,
    get_gui_metadata, get_parent, get_systems # MTK internals
using PlantGraphs, MultiScaleTreeGraph # graphs
using RecipesBase # visualisation

# # Include src files
include("graph_nodetypes.jl")
include("graph_functions.jl")
include("graph_reading.jl")
include("graph_conversion.jl")
include("generate_system.jl")
include("shapes.jl")
include("smoothfuncs.jl")
include("func_modules.jl")
include("plotting.jl")

end # module