module PlantModules

export generate_system, alter_defaults, plotgraph, plotnode, readXEG, convert_to_MTG, convert_to_PG

using ModelingToolkit, DifferentialEquations, Unitful # Simulation packages
import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv
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