@test PlantModules.plotnode(sol, PlantModules.getnodes(plantstructure)[1]) isa Vector{AbstractPlot}
@test PlantModules.plotnode(sol, PlantModules.getnodes(plantstructure)[5], varname = :N) isa AbstractPlot

@test PlantModules.plotgraph(sol, plantstructure) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, plantstructure, structmod = :Forest) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, plantstructure, varname = :ΣF_P) isa AbstractPlot
@test PlantModules.plotgraph(sol, plantstructure, structmod = :Grassland, varname = :P) isa AbstractPlot

@test PlantModules.plotgraph(sol, plantstructure) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, plantstructure, structmod = :Forest) isa Vector{AbstractPlot}
@test PlantModules.plotgraph(sol, plantstructure, varname = :ΣF_P) isa AbstractPlot
@test PlantModules.plotgraph(sol, plantstructure, structmod = :Grassland, varname = :P) isa AbstractPlot
