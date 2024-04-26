if split(pwd(), '\\')[end] != "tutorials"
    cd("./tutorials")
end
using Pkg; Pkg.activate(".")
include("../src/PlantModules.jl"); using .PlantModules
using PlantGraphs, ModelingToolkit, DifferentialEquations, Unitful, Plots, MultiScaleTreeGraph
import GLMakie.draw

# Pkg.develop(path = "..")

using Pluto
Pluto.run()