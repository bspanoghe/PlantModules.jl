module PlantModules

export generate_system, plotgraph, plotnode

using ModelingToolkit
import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv
using DifferentialEquations
using Unitful
using PlantGraphs
using Plots

include("graphfuncs.jl")
include("generate_system.jl")
include("shapes.jl")
include("func_modules.jl")
include("plotting.jl")

end # module