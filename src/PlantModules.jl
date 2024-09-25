module PlantModules

export Shape, Sphere, Cilinder, Cuboid
export volume, cross_area, surface_area
export generate_system, PlantStructure, PlantFunctionality, alter_defaults
export hydraulic_module, constant_carbon_module, environmental_module, Ψ_soil_module, Ψ_air_module, sizedep_K_module, constant_K_module, hydraulic_connection
export readXEG, convert_to_MTG, convert_to_PG
export plotgraph, plotnode

using ModelingToolkit, DifferentialEquations, Unitful # Simulation packages
import ModelingToolkit: get_eqs, get_unknowns, get_ps, get_parameter_dependencies, get_observed, get_continuous_events, get_discrete_events, get_defaults, get_systems, get_name, get_iv, get_gui_metadata
using PlantGraphs, MultiScaleTreeGraph # Graph packages
using Plots # Visualisation packages

include("graph_nodetypes.jl")
include("graph_functions.jl")
include("graph_reading.jl")
include("graph_conversion.jl")
include("generate_system.jl")
include("shapes.jl")
include("func_modules.jl")
include("plotting.jl")
include("qualityoflife.jl")

end # module