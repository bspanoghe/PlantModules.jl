# System creation

```@docs
generate_system
PlantStructure
PlantCoupling
PlantParameters
```

# Shapes

```@docs
ModuleShape
Sphere
Cylinder
Cuboid
```

# Graph functions
```@docs
getnodes
getneighbors
getattributes
getstructmod
getid

getroot
getvariables
getchildren
getparent
```

# Shape functions
```@docs
getdimensionality
correctdimensionality
volume
cross_area
surface_area
```

# Smooth functions
```@docs
logsumexp
smooth_daynight
```

# Node modules
```@docs
hydraulic_module
environmental_module
constant_carbon_module
simple_photosynthesis_module
Ψ_soil_module
Ψ_air_module
K_module
constant_K_module
```

# Edge modules
```@docs
hydraulic_connection
constant_hydraulic_connection
daynight_hydraulic_connection
```

# Graph reading and converting
```@docs
readXEG
convert_to_MTG
convert_to_PG
```

# System remaking
```@docs
remake_graphsystem
remake_graphsystem!
get_subsystem_variables
```

# Plotting
```@docs
plotstructure
plotgraph
plotnode
```