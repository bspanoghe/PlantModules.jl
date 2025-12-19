testgraph = readXEG("./testdata/test.xeg")

@test testgraph isa Dict
@test PlantModules.getnodes(testgraph) isa Vector{Dict}

coolnode_id = 5069719
coolnode = testgraph[coolnode_id]

@test coolnode isa Dict{Symbol, Any}
@test PlantModules.getid(coolnode) == coolnode_id

@test PlantModules.getattributes(coolnode) == Dict(:angle => -19.989698)
@test PlantModules.getstructmod(coolnode) == :RU

coolnode_nbs = PlantModules.getneighbors(coolnode, testgraph)
@test issetequal(PlantModules.getid.(coolnode_nbs), [5069718, 5069720, 5069722])
