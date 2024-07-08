testgraph = readXEG("./testdata/test.xeg")

@test testgraph isa Dict
@test PlantModules.nodes(testgraph) isa Vector{Dict}

coolnode_id = 5069719
coolnode = testgraph[coolnode_id]

@test coolnode isa Dict{Symbol, Any}
@test PlantModules.id(coolnode) == coolnode_id

@test PlantModules.attributes(coolnode) == Dict(:angle => -19.989698)
@test PlantModules.structmod(coolnode) == :RU

coolnode_nbs = PlantModules.neighbours(coolnode, testgraph)
@test issetequal(PlantModules.id.(coolnode_nbs), [5069718, 5069720, 5069722])