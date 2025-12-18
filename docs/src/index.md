# PlantModules.jl

`PlantModules.jl` is a modular modelling framework built on top of `ModelingToolkit.jl` created to simulate functional-structural plant models. Important features of our package include:
1. Extensibility: the entire package is designed to be easily extended and customized.
1. Speed: simulations are fast by virtue of being embedded in the SciML modeling ecosystem.
1. Base functionality: the package includes a base layer of functionality based on water relations that allows it to work out of the box.
1. Clear model representation: Differential-algebraic systems of equations provide an intuitive representation of functional processes.

## Getting started

To use `PlantModules.jl`, install the package as follows:
```julia
using Pkg; Pkg.add("PlantModules.jl")
```
To learn how to use the package, please visit the [Examples]("./tutorials").

Please note that the package is in an early version and the documentation is still under active construction.