@test PlantModules.plotgraph(sol, plantstructure, PlantModules.getnodes(plantstructure)[1:3], varname = :N) isa AbstractPlot
@test PlantModules.plotgraph(sol, plantstructure, varname = :ΣF_P) isa AbstractPlot
@test PlantModules.plotgraph(sol, plantstructure, structmod = :Grassland, varname = :ΣF_P) isa AbstractPlot