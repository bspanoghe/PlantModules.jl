module PlantModules

export Shape, Sphere, Cylinder, Cuboid
export volume, cross_area, surface_area
export generate_system, PlantStructure, PlantFunctionality
export hydraulic_module, constant_carbon_module, environmental_module, Ψ_soil_module, Ψ_air_module
export K_module, constant_K_module, hydraulic_connection, const_hydraulic_connection, evaporation_connection
export readXEG, convert_to_MTG, convert_to_PG
export plotgraph, plotnode

using SciMLBase # Interface with ODE solver packages
using ModelingToolkit, Unitful # Simulation
import ModelingToolkit: get_eqs, get_unknowns, get_ps, get_parameter_dependencies, get_observed, get_continuous_events, get_discrete_events, get_defaults, get_systems, get_name, get_iv, get_gui_metadata
using PlantGraphs, MultiScaleTreeGraph # Graphs
using RecipesBase # Visualisation

include("graph_nodetypes.jl")
include("graph_functions.jl")
include("graph_reading.jl")
include("graph_conversion.jl")
include("generate_system.jl")
include("shapes.jl")
include("func_modules.jl")
include("plotting.jl")

end # module