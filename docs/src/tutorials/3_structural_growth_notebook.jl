### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ efdd1bf3-ee36-447d-9741-eee35583c125
using Pkg

# ╔═╡ 98cb54a9-bfa8-4f9a-bdd3-f211a98c975b
begin
	Pkg.activate()
	using Revise
end

# ╔═╡ 399e95e9-be88-4a0b-8c88-03d74fe0d76b
begin
	Pkg.activate("../..")
	using PlantModules
	using ModelingToolkit, OrdinaryDiffEq, Plots
	using VirtualPlantLab, GLMakie
	import VirtualPlantLab: Mesh
end

# ╔═╡ 2dfdb97f-361c-47bd-b65a-41b02bc8bc57
md"# Tutorial 3: Functional-structural growth modelling"

# ╔═╡ 87a39d4c-6f22-447e-bd45-1f133b5bb7b3
md"""
In the tutorials so far, we have only considered the simulation of water flows for static plant structures. In this tutorial, we will finally create a real FSPM by simulating functional- and structural growth simultaneously, using [`VirtualPlantLab.jl`](https://virtualplantlab.com/stable/) for the simulation of structural growth.
"""

# ╔═╡ e52a8761-f703-4dec-b02b-62d3cc831b4f
md"## Setup"

# ╔═╡ a36a2a5b-de0b-46a5-a97c-d0c6b1f2f9f4
md"## Context"

# ╔═╡ 25160d35-aa6b-4ed5-9199-a54d882dbfc9
md"""
We build further on the [tree growth modelling tutorial](https://virtualplantlab.com/stable/tutorials/from_tree_forest/tree/) from `VirtualPlantLab.jl`. 

In the original version, the tree grows by two rewriting steps:
- Tree meristems grow into phytomers, occuring every rewrite step,
- Buds grow branches, occuring with a probability proportional to the amount of phytomers to the apical meristem,
and one elongation step:
- Every internode elongates for a set fraction of its current length.

In this tutorial, we will replace the basic exponential elongation with growth based on water relations, and make the structural changes depend on the plant's hydraulic status.
"""

# ╔═╡ 441f1768-b9b3-49e1-87ab-060ff720948f
md"## Original tutorial"

# ╔═╡ 83fadc9e-74a5-4ddd-915a-0b8692c9280f
md"### Structural module definition"

# ╔═╡ b27f3786-84b0-469f-ab0c-a60f02030efd
Base.@kwdef struct Meristem <: VirtualPlantLab.Node 
	D = [0.5] # 
end

# ╔═╡ b99eb3ee-26dc-4396-9a61-aba45f7462cf
Base.@kwdef struct Bud <: VirtualPlantLab.Node
	D = [0.5]
end

# ╔═╡ 43bb6044-e567-4f20-8be0-5108f5deead8
Base.@kwdef struct Node <: VirtualPlantLab.Node
	D = [0.5]
end

# ╔═╡ 4ea2dfcf-cc13-46f1-8f70-f42af7a43a09
Base.@kwdef struct BudNode <: VirtualPlantLab.Node
	D = [0.5]
end

# ╔═╡ fc109086-c42c-474d-8ab0-4add40a07eed
Base.@kwdef mutable struct Internode <: VirtualPlantLab.Node
	D = [0.5, 10] # internodes are cylinders starting at 0.5 cm radius, 10 cm length
end

# ╔═╡ 2ff96078-4a91-427b-8c39-689d47600c05
Base.@kwdef struct Leaf <: VirtualPlantLab.Node
	D = [20, 10, 0.05] # leaves are cuboids starting at 20 cm long, 10 cm wide and 0.05 cm thick
end

# ╔═╡ fea00274-c739-48cb-8085-4f5713ebeb9c
md"### Parameter definition"

# ╔═╡ 120217e8-c708-4a49-8f2c-1191690ce0db
Base.@kwdef struct treeparams
	growth::Float64 = 0.1
	budbreak::Float64 = 0.25
	phyllotaxis::Float64 = 140.0
	leaf_angle::Float64 = 30.0
	branch_angle::Float64 = 45.0
end

# ╔═╡ 5a1e5b8a-be42-45e8-9abc-5cecd9dc83ef
md"### Visualisation functions"

# ╔═╡ 1f2b0dc8-84db-4612-841a-6ddaf08345e5
function VirtualPlantLab.feed!(turtle::Turtle, i::Internode, vars)
    # Rotate turtle around the head to implement elliptical phyllotaxis
    rh!(turtle, vars.phyllotaxis)
    HollowCylinder!(
		turtle, length = i.D[2], height = i.D[1],
		width = i.D[1], move = true, colors = RGB(0.5,0.4,0.0)
	)
    return nothing
end

# ╔═╡ e062618b-3913-4dd8-8ef4-fdcba785fa61
function VirtualPlantLab.feed!(turtle::Turtle, l::Leaf, vars)
    ra!(turtle, -vars.leaf_angle)
    Ellipse!(turtle, length = l.D[1], width = l.D[2], move = false,
			 colors = RGB(0.2, 0.6, 0.2))
    ra!(turtle, vars.leaf_angle)
	
    return nothing
end

# ╔═╡ d6fa9897-05ae-4388-a299-c08cb62c1414
function VirtualPlantLab.feed!(turtle::Turtle, b::BudNode, vars)
    ra!(turtle, -vars.branch_angle)
end

# ╔═╡ 2c04e08b-1006-4e88-b754-1297eed4dd9c
md"### Structural growth rules"

# ╔═╡ eb87190d-26e5-4e37-bb09-15b1e89da868
meristem_rule = Rule(
	Meristem,
	rhs = mer -> Node() + (Bud(), Leaf()) + Internode() + Meristem()
)

# ╔═╡ d0e32e1a-733e-4791-82ec-c4a20745a351
function prob_break(bud)
    node = parent(bud)
    check, steps = has_descendant(
		node, 
		condition = n -> data(n) isa Meristem
	)
    steps = Int(ceil(steps/2))
    if check
        prob =  min(1.0, steps*graph_data(bud).budbreak)
        return rand() < prob
    else
        error("No meristem found in branch")
    end
end

# ╔═╡ fb31e01e-efda-4e72-828c-efd5a286d673
branch_rule = Rule(
	Bud,
	lhs = prob_break,
	rhs = bud -> BudNode() + Internode() + Meristem()
)

# ╔═╡ 421e7e84-3c51-4e64-adc4-2a98a4536844
axiom = Internode() + Meristem()

# ╔═╡ 04d3fc8a-2ab7-4085-8fd6-4ab3f09eeebe
tree = Graph(axiom = axiom, rules = (meristem_rule, branch_rule), data = treeparams())

# ╔═╡ 1da39117-7f6d-4907-aad7-bcc0d0fc4be4
getInternode = Query(Internode)

# ╔═╡ 1bb9f437-44fb-4a7c-9eee-2730d1487e87
function elongate!(tree, query)
    for x in apply(tree, query)
        x.D = x.D .* [1.0 , 1.0 + data(tree).growth]
    end
end

# ╔═╡ 0fda99f0-ac5e-4080-bc44-3378ef41c856
function growth!(tree, query)
    elongate!(tree, query)
    rewrite!(tree)
end

# ╔═╡ 7e7fdcd0-e59e-471a-817a-87a8403440e4
function simulate(tree, query, nsteps)
    new_tree = deepcopy(tree)
    for i in 1:nsteps
        growth!(new_tree, query)
    end
    return new_tree
end

# ╔═╡ e6bdd4a1-135d-4c1b-be8e-24d5e1a3f3b1
newtree = simulate(tree, getInternode, 5)

# ╔═╡ e0ca99d2-84e1-47e3-93bc-324465f5abd3
render(Mesh(newtree))

# ╔═╡ bf56547c-3d1d-4134-ba92-5f6c521e1d07
md"## Extending the model"

# ╔═╡ 7edd2626-a111-45bb-a62c-a526d8698f86
md"### Defining the environment"

# ╔═╡ 1ced3c28-55be-4f4e-99d0-c38f6cc79396
struct Soil <: VirtualPlantLab.Node end

# ╔═╡ 082576ba-0932-432b-b173-c844fb57bfc0
struct Air <: VirtualPlantLab.Node end

# ╔═╡ 3b249efd-55cb-4a65-a245-07296809c6b6
graphs = [newtree, Soil(), Air()];

# ╔═╡ 6461c5f0-5860-433d-8627-53a14dd342e4
intergraph_connections = [(1, 2) => (getnodes(newtree)[1], :Soil), (1, 3) => (:Leaf, :Air)];

# ╔═╡ 52aa8d2e-16f5-49fc-86fd-0f1c3a325b3d
plantstructure = PlantStructure(graphs, intergraph_connections);

# ╔═╡ ca063cc1-cd0e-48de-b65b-0e42c0f0df65
plotstructure(plantstructure)

# ╔═╡ ac0e93c0-2662-4ea4-bcc3-3067aec41d28
md"### Functional definition"

# ╔═╡ 62257b82-ac54-4360-951c-9192edb619a3
md"#### Coupling"

# ╔═╡ 767cb212-1e7d-4909-8be2-d873962c6294
import PlantModules: t, d

# ╔═╡ 9232633a-5df0-4825-b6a9-d7efd5e036d7
function inactive_hydraulic_module(; name)
	@variables ΣF(t)
	eqs = [d(ΣF) ~ 0]

	return System(eqs, t; name)
end

# ╔═╡ 396c9109-d930-4f89-b280-ab23b918391f
function inactive_hydraulic_connection(; name)
    @variables F(t) = 0.0

    eqs = [d(F) ~ 0]
    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = []
	
    return System(eqs, t; name), get_connection_eqset
end

# ╔═╡ 251fd6eb-627b-4560-b692-d02bef8f089c
module_coupling = Dict(
	:Meristem => [inactive_hydraulic_module],
	:Bud => [hydraulic_module, constant_carbon_module, constant_K_module],
    :Node => [hydraulic_module, constant_carbon_module, constant_K_module],
	:BudNode => [hydraulic_module, constant_carbon_module, constant_K_module],
	:Internode => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
);

# ╔═╡ 812ec2a4-b3ec-4240-a7bd-b80ddea7a748
connecting_modules = Dict(
	(:Soil, :Internode) => constant_hydraulic_connection,
	(:Internode, :Node) => hydraulic_connection,
	(:Node, :Bud) => hydraulic_connection,
	(:Node, :BudNode) => hydraulic_connection,
	(:BudNode, :Internode) => hydraulic_connection,
	(:Node, :Leaf) => hydraulic_connection,
	(:Internode, :Meristem) => inactive_hydraulic_connection,
	(:Leaf, :Air) => daynight_hydraulic_connection
);

# ╔═╡ 1f90f558-272b-41b6-a4ba-c7bbcea3cf92
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ 8f986516-75a7-46b2-86f5-4836d8b54193
md"#### Parameters"

# ╔═╡ 6452a19d-fe5a-48b6-8ec8-d8759407a4dc
default_changes = Dict([:Ψ => PlantModules.soilfunc(0.8)])

# ╔═╡ 1d8aba86-aba8-4567-a3f6-648d3b74e8e3
module_defaults = Dict(
	:Bud => Dict(:shape => PlantModules.Sphere()),
	:Node => Dict(:shape => PlantModules.Sphere()),
	:BudNode => Dict(:shape => PlantModules.Sphere()),
	:Leaf => Dict(:shape => PlantModules.Cuboid()),
	:Soil => Dict(:W_max => 1e4, :K => 1.0),
	:Air => Dict(:W_r => 0.7, :K => 1e-3)
);

# ╔═╡ 2f54285c-50cf-49e4-b2b0-cf9aa5fdf585
plantparams = PlantParameters(; default_changes, module_defaults);

# ╔═╡ 25c85fd1-1921-4707-8185-735751c0d914
md"### Creating the system"

# ╔═╡ 8f635abb-c9ea-45ca-81ae-e679e0ba923d
system = generate_system(plantstructure, plantcoupling, plantparams);

# ╔═╡ 7b8bbd0f-c650-49ef-b8d9-c42cd6d8da9e
prob = ODEProblem(system, [], (0.0, 480.0), sparse = true)

# ╔═╡ 8363c723-74e0-49bf-88ae-e164919f4681
sol = solve(prob, FBDF())

# ╔═╡ 64724bb9-05c2-448c-8f72-edb88edbce45
plotgraph(sol, plantstructure, varname = :V, structmod = :Internode)

# ╔═╡ 8743b00d-5699-4201-9a23-0061552495a0
plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Internode, :Bud, :Soil])

# ╔═╡ e48e51dc-de57-4186-891f-7c74fae1ebd1
plotgraph(sol, plantstructure, varname = :W, structmod = :Bud)

# ╔═╡ 853bbb05-d425-4ba4-aa23-c0af0ed4d021
plotgraph(sol, plantstructure, varname = :W, structmod = :Soil)

# ╔═╡ c320aed9-7086-4dcc-8030-b3f281f5a1ec
md"## Speed benchmarking"

# ╔═╡ Cell order:
# ╟─2dfdb97f-361c-47bd-b65a-41b02bc8bc57
# ╟─87a39d4c-6f22-447e-bd45-1f133b5bb7b3
# ╟─e52a8761-f703-4dec-b02b-62d3cc831b4f
# ╠═efdd1bf3-ee36-447d-9741-eee35583c125
# ╠═98cb54a9-bfa8-4f9a-bdd3-f211a98c975b
# ╠═399e95e9-be88-4a0b-8c88-03d74fe0d76b
# ╟─a36a2a5b-de0b-46a5-a97c-d0c6b1f2f9f4
# ╟─25160d35-aa6b-4ed5-9199-a54d882dbfc9
# ╟─441f1768-b9b3-49e1-87ab-060ff720948f
# ╟─83fadc9e-74a5-4ddd-915a-0b8692c9280f
# ╠═b27f3786-84b0-469f-ab0c-a60f02030efd
# ╠═b99eb3ee-26dc-4396-9a61-aba45f7462cf
# ╠═43bb6044-e567-4f20-8be0-5108f5deead8
# ╠═4ea2dfcf-cc13-46f1-8f70-f42af7a43a09
# ╠═fc109086-c42c-474d-8ab0-4add40a07eed
# ╠═2ff96078-4a91-427b-8c39-689d47600c05
# ╟─fea00274-c739-48cb-8085-4f5713ebeb9c
# ╠═120217e8-c708-4a49-8f2c-1191690ce0db
# ╟─5a1e5b8a-be42-45e8-9abc-5cecd9dc83ef
# ╠═1f2b0dc8-84db-4612-841a-6ddaf08345e5
# ╠═e062618b-3913-4dd8-8ef4-fdcba785fa61
# ╠═d6fa9897-05ae-4388-a299-c08cb62c1414
# ╟─2c04e08b-1006-4e88-b754-1297eed4dd9c
# ╠═eb87190d-26e5-4e37-bb09-15b1e89da868
# ╠═d0e32e1a-733e-4791-82ec-c4a20745a351
# ╠═fb31e01e-efda-4e72-828c-efd5a286d673
# ╠═421e7e84-3c51-4e64-adc4-2a98a4536844
# ╠═04d3fc8a-2ab7-4085-8fd6-4ab3f09eeebe
# ╠═1da39117-7f6d-4907-aad7-bcc0d0fc4be4
# ╠═1bb9f437-44fb-4a7c-9eee-2730d1487e87
# ╠═0fda99f0-ac5e-4080-bc44-3378ef41c856
# ╠═7e7fdcd0-e59e-471a-817a-87a8403440e4
# ╠═e6bdd4a1-135d-4c1b-be8e-24d5e1a3f3b1
# ╠═e0ca99d2-84e1-47e3-93bc-324465f5abd3
# ╟─bf56547c-3d1d-4134-ba92-5f6c521e1d07
# ╟─7edd2626-a111-45bb-a62c-a526d8698f86
# ╠═1ced3c28-55be-4f4e-99d0-c38f6cc79396
# ╠═082576ba-0932-432b-b173-c844fb57bfc0
# ╠═3b249efd-55cb-4a65-a245-07296809c6b6
# ╠═6461c5f0-5860-433d-8627-53a14dd342e4
# ╠═52aa8d2e-16f5-49fc-86fd-0f1c3a325b3d
# ╠═ca063cc1-cd0e-48de-b65b-0e42c0f0df65
# ╟─ac0e93c0-2662-4ea4-bcc3-3067aec41d28
# ╟─62257b82-ac54-4360-951c-9192edb619a3
# ╠═767cb212-1e7d-4909-8be2-d873962c6294
# ╠═9232633a-5df0-4825-b6a9-d7efd5e036d7
# ╠═396c9109-d930-4f89-b280-ab23b918391f
# ╠═251fd6eb-627b-4560-b692-d02bef8f089c
# ╠═812ec2a4-b3ec-4240-a7bd-b80ddea7a748
# ╠═1f90f558-272b-41b6-a4ba-c7bbcea3cf92
# ╟─8f986516-75a7-46b2-86f5-4836d8b54193
# ╠═6452a19d-fe5a-48b6-8ec8-d8759407a4dc
# ╠═1d8aba86-aba8-4567-a3f6-648d3b74e8e3
# ╠═2f54285c-50cf-49e4-b2b0-cf9aa5fdf585
# ╟─25c85fd1-1921-4707-8185-735751c0d914
# ╠═8f635abb-c9ea-45ca-81ae-e679e0ba923d
# ╠═7b8bbd0f-c650-49ef-b8d9-c42cd6d8da9e
# ╠═8363c723-74e0-49bf-88ae-e164919f4681
# ╠═64724bb9-05c2-448c-8f72-edb88edbce45
# ╠═8743b00d-5699-4201-9a23-0061552495a0
# ╠═e48e51dc-de57-4186-891f-7c74fae1ebd1
# ╠═853bbb05-d425-4ba4-aa23-c0af0ed4d021
# ╟─c320aed9-7086-4dcc-8030-b3f281f5a1ec
