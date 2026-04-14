### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ddb74042-dbe4-4a1d-94e7-81f4af878af6
using Pkg; Pkg.activate("../..")

# ╔═╡ 76faf1ec-1c1e-46af-929c-ad0ab80485c6
using PlantModules

# ╔═╡ 5f1e0631-1df3-4ae4-8650-736556ad7f10
using ModelingToolkit, OrdinaryDiffEq, Plots

# ╔═╡ 368f9b09-1c40-421e-9557-9641d0b3bf5c
using VirtualPlantLab, ColorTypes, GLMakie

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
struct Meristem <: VirtualPlantLab.Node end

# ╔═╡ b99eb3ee-26dc-4396-9a61-aba45f7462cf
struct Bud <: VirtualPlantLab.Node end

# ╔═╡ 43bb6044-e567-4f20-8be0-5108f5deead8
struct Node <: VirtualPlantLab.Node end

# ╔═╡ 4ea2dfcf-cc13-46f1-8f70-f42af7a43a09
struct BudNode <: VirtualPlantLab.Node end

# ╔═╡ fc109086-c42c-474d-8ab0-4add40a07eed
Base.@kwdef mutable struct Internode <: VirtualPlantLab.Node
	length::Float64 = 0.10 # Internodes start at 10 cm
end

# ╔═╡ 2ff96078-4a91-427b-8c39-689d47600c05
Base.@kwdef struct Leaf <: VirtualPlantLab.Node
	length::Float64 = 0.20 # Leaves are 20 cm long
	width::Float64  = 0.1 # Leaves are 10 cm wide
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
    HollowCylinder!(turtle, length = i.length, height = i.length/15, width = i.length/15,
                move = true, colors = RGB(0.5,0.4,0.0))
    return nothing
end

# ╔═╡ e062618b-3913-4dd8-8ef4-fdcba785fa61
# Create geometry + color for the leaves
function VirtualPlantLab.feed!(turtle::Turtle, l::Leaf, vars)
    # Rotate turtle around the arm for insertion angle
    ra!(turtle, -vars.leaf_angle)
    # Generate the leaf
    Ellipse!(turtle, length = l.length, width = l.width, move = false,
             colors = RGB(0.2, 0.6, 0.2))
    # Rotate turtle back to original direction
    ra!(turtle, vars.leaf_angle)
    return nothing
end

# ╔═╡ d6fa9897-05ae-4388-a299-c08cb62c1414
# Insertion angle for the bud nodes
function VirtualPlantLab.feed!(turtle::Turtle, b::BudNode, vars)
    # Rotate turtle around the arm for insertion angle
    ra!(turtle, -vars.branch_angle)
end

# ╔═╡ 2c04e08b-1006-4e88-b754-1297eed4dd9c
md"### Structural growth rules"

# ╔═╡ eb87190d-26e5-4e37-bb09-15b1e89da868
meristem_rule = Rule(
	TreeTypes.Meristem,
	rhs = mer -> TreeTypes.Node() + 
		(TreeTypes.Bud(), TreeTypes.Leaf()) +
		TreeTypes.Internode() + TreeTypes.Meristem()
)

# ╔═╡ d0e32e1a-733e-4791-82ec-c4a20745a351
function prob_break(bud)
    # We move to parent node in the branch where the bud was created
    node =  parent(bud)
    # We count the number of internodes between node and the first Meristem
    # moving down the graph
    check, steps = has_descendant(node, condition = n -> data(n) isa TreeTypes.Meristem)
    steps = Int(ceil(steps/2)) # Because it will count both the nodes and the internodes
    # Compute probability of bud break and determine whether it happens
    if check
        prob =  min(1.0, steps*graph_data(bud).budbreak)
        return rand() < prob
    # If there is no meristem, an error happened since the model does not allow for this
    else
        error("No meristem found in branch")
    end
end

# ╔═╡ fb31e01e-efda-4e72-828c-efd5a286d673
branch_rule = Rule(
	TreeTypes.Bud,
	lhs = prob_break,
	rhs = bud -> TreeTypes.BudNode() +
		TreeTypes.Internode() + TreeTypes.Meristem()
)

# ╔═╡ 421e7e84-3c51-4e64-adc4-2a98a4536844
axiom = TreeTypes.Internode() + TreeTypes.Meristem()

# ╔═╡ 04d3fc8a-2ab7-4085-8fd6-4ab3f09eeebe
tree = Graph(axiom = axiom, rules = (meristem_rule, branch_rule), data = TreeTypes.treeparams())

# ╔═╡ 1da39117-7f6d-4907-aad7-bcc0d0fc4be4
getInternode = Query(TreeTypes.Internode)

# ╔═╡ 1bb9f437-44fb-4a7c-9eee-2730d1487e87
function elongate!(tree, query)
    for x in apply(tree, query)
        x.length = x.length*(1.0 + data(tree).growth)
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
newtree = simulate(tree, getInternode, 15)

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
intergraph_connections = [(1, 2) => (getnodes(plant)[1], :Soil), (1, 3) => (:BranchTip, :Air)];

# ╔═╡ 52aa8d2e-16f5-49fc-86fd-0f1c3a325b3d
plantstructure = PlantStructure(graphs, intergraph_connections);

# ╔═╡ ca063cc1-cd0e-48de-b65b-0e42c0f0df65
plotstructure(plantstructure, names = "")

# ╔═╡ ac0e93c0-2662-4ea4-bcc3-3067aec41d28
md"### Functional definition"

# ╔═╡ 251fd6eb-627b-4560-b692-d02bef8f089c
module_coupling = Dict(
	:Meristem => [hydraulic_module, constant_carbon_module, K_module],
	:Bud => [],
    :Node => [hydraulic_module, constant_carbon_module, K_module],
	:BudNode => [hydraulic_module, needle_area_module,
				   constant_carbon_module, K_module],
	:Internode => [],
	:Leaf => [],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module],
);

# ╔═╡ 812ec2a4-b3ec-4240-a7bd-b80ddea7a748
connecting_modules = Dict(
	(:Soil, :Stem) => constant_hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Branch) => hydraulic_connection,
    (:Branch, :Branch) => hydraulic_connection,
	(:Branch, :BranchTip) => hydraulic_connection,
	(:BranchTip, :Air) => fixed_transpiration_connection,
);

# ╔═╡ 1f90f558-272b-41b6-a4ba-c7bbcea3cf92
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ c320aed9-7086-4dcc-8030-b3f281f5a1ec
md"## Speed benchmarking"

# ╔═╡ Cell order:
# ╟─2dfdb97f-361c-47bd-b65a-41b02bc8bc57
# ╟─87a39d4c-6f22-447e-bd45-1f133b5bb7b3
# ╟─e52a8761-f703-4dec-b02b-62d3cc831b4f
# ╠═ddb74042-dbe4-4a1d-94e7-81f4af878af6
# ╠═76faf1ec-1c1e-46af-929c-ad0ab80485c6
# ╠═5f1e0631-1df3-4ae4-8650-736556ad7f10
# ╠═368f9b09-1c40-421e-9557-9641d0b3bf5c
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
# ╠═251fd6eb-627b-4560-b692-d02bef8f089c
# ╠═812ec2a4-b3ec-4240-a7bd-b80ddea7a748
# ╠═1f90f558-272b-41b6-a4ba-c7bbcea3cf92
# ╟─c320aed9-7086-4dcc-8030-b3f281f5a1ec
