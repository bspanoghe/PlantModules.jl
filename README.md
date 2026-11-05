[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://bspanoghe.github.io/PlantModules.jl/stable)

# PlantModules.jl
`PlantModules.jl` is a functional-structural plant modeling framework implemented in Julia, created on top of `ModelingToolkit.jl`. Some highlights:
1. Extensibility: the entire package is designed to be easily extended and customized.
2. Speed: simulations are fast by virtue of being embedded in the SciML modeling ecosystem.
3. Base functionality: the package includes a base set of functional processes describing hydraulics and hydraulics-driven growth, allowing the creation of plant models out of the box.
4. Clear model representation: Differential-algebraic systems of equations provide an intuitive representation of functional processes.
5. All Julia: As a pure Julia implementation, it integrates seamlessly with many other Julia packages such as the plant modeling framework [VirtualPlantLab.jl](https://github.com/VirtualPlantLab/VirtualPlantLab.jl).

## Installation
The package can be installed as follows, in the Julia REPL:
```julia
] add PlantModules
```

## Usage
Please see [the documentation page](https://bspanoghe.github.io/PlantModules.jl/stable) for an explanation on how to use the package.

## Citation
If you use this package in a publication, please cite the following paper:
[PlantModules. jl: a framework for modular plant modelling](https://doi.org/10.1093/insilicoplants/diag016)

BibTeX:
```
@article{spanoghe2026,
	title = {{PlantModules}.jl: a framework for modular plant modelling},
	shorttitle = {{PlantModules}.jl},
	author = {Spanoghe, Bram and De Swaef, Tom and De Baets, Bernard and Stock, Michiel},
	journal = {in silico Plants},
	volume = {8},
	number = {2},
	pages = {diag016},
	publisher = {Oxford Academic},
	month = jul,
	year = {2026},
	url = {https://doi.org/10.1093/insilicoplants/diag016},
	doi = {10.1093/insilicoplants/diag016},
}
```
