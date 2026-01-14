# PlantModules.jl
`PlantModules.jl` is a functional-structural plant modeling framework implemented in Julia, created on top of `ModelingToolkit.jl`. Some highlights:
1. Extensibility: the entire package is designed to be easily extended and customized.
2. Speed: simulations are fast by virtue of being embedded in the SciML modeling ecosystem.
3. Base functionality: the package includes a base set of functional processes describing hydraulics and hydraulics-driven growth, allowing the creation of plant models out of the box.
4. Clear model representation: Differential-algebraic systems of equations provide an intuitive representation of functional processes.
5. All Julia: As a pure Julia implementation, it integrates seamlessly with many other Julia packages such as the plant modeling framework [VirtualPlantLab.jl](https://github.com/VirtualPlantLab/VirtualPlantLab.jl).

## Installation
The package is currently not yet in the Julia registry. Therefore, installation of the package can be done as follows, in the Julia REPL:
```julia
] add https://github.com/bspanoghe/PlantModules.jl.git
```

## Usage
To learn how to use the package, please visit the [Examples]("./tutorials").

Please note that the package is in an early version and the documentation is still under active construction.