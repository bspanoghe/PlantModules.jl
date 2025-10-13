### A Pluto.jl notebook ###
# v0.20.14

using Markdown
using InteractiveUtils

# ╔═╡ 3fc5da40-acde-4398-8a53-99794e42ba94
using Pkg; Pkg.activate("../..")

# ╔═╡ 04361839-2e2f-4475-8830-277c4b07e3f6
using PlantModules

# ╔═╡ 45991fe2-1dae-4633-b3c2-81b3e076c1eb
using PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Unitful, Plots

# ╔═╡ 8adf74c7-7fe0-42fd-bde5-9942f30fea36
using DataInterpolations, Measurements

# ╔═╡ df79e96c-a021-4a24-ae8c-4090d014c1f8
using PlutoUI; TableOfContents()

# ╔═╡ e04d4d44-3795-49d0-91d3-645ff3c8265e
md"# Tutorial 2: Package validation"

# ╔═╡ 61d4d897-e1ad-4f04-b0be-1df409c4e69f
md"""
In this tutorial, we cover some of the more advanced functionality from the package in order to recreate an hydraulic tree growth model from literature. The following topics are discussed:
- Defining a new geometrical shape for the structural modules.
- Defining a large plant structure with the use of graph rewriting (which can be considered a generalisation of L-systems).
- Defining new functional modules.
    - Using time-series data in a functional module.
    - Defining a new, asymmetric connection module.
- Using the advanced visualisation functionality to compare simulation outputs with data.
"""

# ╔═╡ 21271cd8-eef1-485a-978d-f46afdb753b6
md"## Setup"

# ╔═╡ b2f4ec18-5447-4999-9da6-52ca770a8120
md"## Problem definition"

# ╔═╡ ac068794-a6cc-4563-899a-4b5ca1b2c22b
md"""
In this tutorial, we recreate [the dynamic sap flow model by Perämäki et al. (2001)](https://doi.org/10.1093/treephys/21.12-13.889) in `PlantModules.jl`. The model simulates stem diameter variation in Scots pine trees caused by transpiration-driven changes in water tension and wood elasticity. The tree is modelled as a series of cylindrical stem and branch segments with a constant length and a radius that decreases as they move farther from the base of the tree. Water flow is modelled as a linear function of the water pressure difference between segments, including the effects of gravity, and the cross-sectional area connecting two segments. For boundary conditions, the outermost branch segments have a set transpiration rate per unit needle area based on measurements, and the base stem segment has a water inflow based on the water potential difference with the soil, of which the water potential is also based on measurements.

This model is an interesting case study for a few reasons:
- It is also based on hydraulics-based plant growth, so we can use the hydraulics-based core functionality from the package.
- The structure is composed of reasonably large number of segments from two different types, which lends itself well to a modular modeling approach.
- There are plenty of interesting things happening both in the structural and functional definition that allow us to showcase the extensibility of the package.
- The paper includes input- and output data, allowing for model validation.
"""

# ╔═╡ 78f2ca6c-9a4e-4b3f-b0f5-d7df2a78fa36
md"## Preparatory calculations"

# ╔═╡ 5a0bf882-68e2-4ad3-ae1f-43cac5d8f890
md"""
In order to replicate the model from the paper, we start by defining the structural and functional parameters used there. The paper was transparent enough that all functional parameters can simply be copied (with some unit conversions), though the structure of the tree was not provided exactly, so some estimations had to be made. The most influential of these on the model output is the fraction of the stem radius which consists of water-conducting sapwood (in contrast to the non-conducting heartwood), which was estimated visually.
"""

# ╔═╡ 1371c204-ba0a-4f18-9df3-d2614277bc11
begin
    # dimensions of compartments

    ## stem segments
    segment_length = 0.2 * 1e2 # cm; from m (input)
    tree_length = 13.2 * 1e2 # cm; from m (input)
    n_segments = tree_length ÷ segment_length

    ## stem radius at tree base
    stem_r0 = 0.2 * 1e2 / 2 # cm; from the diameter in m (input)
    radial_fraction_sapwood = 0.25 # (estimate)
    sapwood_r0 = stem_r0 * radial_fraction_sapwood

    ## crown
    crown_base_height = 6.0 * 1e2 # cm; from m (input)
    crown_length = tree_length - crown_base_height
    n_crown_segments = crown_length ÷ segment_length
    crown_start_age = (n_segments - n_crown_segments)

    ## branching
    branching_frequency = 3 # (estimate)
    n_branches = n_crown_segments ÷ branching_frequency

    Δr = sapwood_r0/(n_branches + 1) 
		# how much the radius of a stem segment decreases whenever the stem branches
		# we assume the tree decreases linearly in radius from the given starting radius at the bottom to a radius of near 0 at the top
		# (+1 to prevent stem radius from hitting 0)

    ## branch radius
    	# calculate using da vinci's rule: total cross area remains constant
		# before and after branching
    function branch_radius(r, Δr)
        area_0 = r^2*pi
        area_new = (r - Δr)^2*pi
        area_branch = area_0 - area_new
        return sqrt(area_branch / pi)
    end

    ## needle areas
    total_needle_area = 50.8 * 1e2^2 # cm^2; from m^2 (input)
    needle_areas = [i * total_needle_area / sum(1:n_branches)
					for i in n_branches:-1:1]

    ## branchgrowth frequency
    branchgrowth_frequency = 3 
		# how many steps it takes branches to grow a new segment (estimate)

    ## root dimensions
    total_root_area = 96.1 * 1e2^2 # cm^2; from m^2 (input)

    # parameters

    ## hydraulic conductivities
    permeability = 2.5e-12 * 1e2^2 # cm^2; from m^2 (input)
    η = 1.0 * 1e-9 * (1/3600) # MPa hr; from mPa s (dynamic viscosity of water)
    l = segment_length # cm
    ρ_w = 1.0 # g / cm^3 (density of water)
    K_s_stem = ρ_w * permeability / η / l # g / hr / MPa / cm^2

    L_p = 3.2e-8 * 1e2 * 3600 # cm / hr / MPa; from m / s / MPa (input)
    K_roots = ρ_w * L_p * total_root_area # g / hr / MPa

    ## elastic modulus
    ϵ_D_stem = [0.15 * 1e3, 17.5 * 0.15 * 1e3] # MPa; from GPa (input)
end;

# ╔═╡ fdedf540-9c04-476f-9ffb-d7bc82e81f79
md"## Structural definition"

# ╔═╡ e26238dd-d323-4161-9e95-bece1af1b3c5
md"""
We introduce two new possibilities for the definition of plant structure. 

Firstly, we need to define a new shape for our structural modules: the stem and branches consist of a hollow cylinder of growing, water-conducting sapwood with an "inert" core of heartwood in the middle. The base hydraulics module models a growing, water-conducting element, so we want to only model the sapwood.

Secondly, the model structure described in the paper consists of a reasonably large amount of segments (over a hundred) that adhere to some basic rules, so we define the structure by use of graph rewriting, using the functionality from [PlantGraphs.jl](https://github.com/VirtualPlantLab/PlantGraphs.jl).
"""

# ╔═╡ 97b3b0b1-0408-46ba-96e0-01eadab25ea0
md"### Defining a new shape"

# ╔═╡ 479acd0d-ce44-411d-97bf-e12a49004a5d
md"""
In order to define a new shape, we define a new composite type as a subtype of the `PlantModules.Shape` abstract type. Afterwards, we extend the following functions for the shape `MyShape` with dimensions the vector `D`, also listed in the `PlantModules.Shape` docstring:
- `getdimensionality(m::MyShape)`: Return how many dimensions define the shape. 
- `volume(m::MyShape, D)`: Calculate the volume of the shape. 
- `cross_area(m::MyShape, D)`: Calculate the cross-sectional area of the shape. 
- `surface_area(m::MyShape, D)`: Calculate the surface area of the shape.
  - This one is used for the photosynthesis module, which we won't be using this tutorial, so we don't need to extend it.
"""

# ╔═╡ 2a4b2b06-a057-4fef-b8fe-30700b39b564
"""
    HollowCylinder <: Shape

A compartment shape representing a hollow cylinder. It is defined by two dimensions: the radius and the length, and the attribute `frac_sapwood`, denoting the fraction of the radius that is conducting sapwood.
"""
struct HollowCylinder <: PlantModules.Shape
    frac_sapwood::Float64
end

# ╔═╡ 6d1a8b68-e8d3-43d9-84a3-8b3fdd600913
PlantModules.getdimensionality(::HollowCylinder) = 2

# ╔═╡ 1ca717b5-7c67-4727-a984-7fd6914a3257
PlantModules.cross_area(hc::HollowCylinder, D::AbstractArray) = D[1]^2 * pi * (2 / hc.frac_sapwood - 1)

# ╔═╡ c0ce423e-ddec-4a2f-b724-8d55d3d05073
PlantModules.volume(hc::HollowCylinder, D::AbstractArray) = cross_area(hc, D) * D[2]

# ╔═╡ 17a090c0-f673-4c2d-b1e4-8f9f8a360179
PlantModules.surface_area(hc::HollowCylinder, D::AbstractArray) = (D[1]/frac_sapwood)^2 * pi * D[2]

# ╔═╡ d59a64b2-9fdf-4723-b764-b27330a050fc
md"### Plant structural module types"

# ╔═╡ c7150b44-966b-4c61-b66d-b8080cd51acd
md"""
The plant structural modules are defined with an attribute for their dimensions and their height, as these are variables that need to differ between nodes. We define separate modules for the outermost segments of both the stem and the branches: this makes the graph rewriting easier, and the tips of the branches have a different functional behaviour from the rest of the branch, as all transpiration is assigned to them. The branch tips have a `needle_area` attribute which describes the total area of needles on that branch, which is used to calculate its transpiration.
"""

# ╔═╡ d7ea529f-4fe2-4daf-9c18-bb4621f0c86d
struct Stem <: PlantGraphs.Node
    D::Vector
    height
end

# ╔═╡ ee63840c-392c-4e4a-88d9-a41ccb88fa6f
struct Branch <: PlantGraphs.Node
    D::Vector
    height
end

# ╔═╡ dc48821b-ffb3-4599-b744-b69364e7ad6f
struct StemTip <: PlantGraphs.Node
    D::Vector
    height
    branch_nr
    age
end

# ╔═╡ 578a128c-37e7-4b4d-8ba1-513565027110
struct BranchTip <: PlantGraphs.Node
    D::Vector
    height
    age
    needle_area
end

# ╔═╡ e319323b-6d26-410b-8574-411686d4067c
md"### Growing the plant"

# ╔═╡ b7d2469c-642e-435d-a399-5406c55e01c4
md"""
Graph rewriting is done by defining a set of rules that dictate how the graph needs to be rewritten every step and an initial graph. The rewriting rules are then applied a number of time to achieve the final graph.

[PlantGraphs.jl](https://github.com/VirtualPlantLab/PlantGraphs.jl) implements their rules with three arguments:
- The node type the rule applies to.
- A function describing the condition for when to apply the rule.
- A function describing how the node is transformed when the rule is applied.

This section is very context-specific and not functionality provided by this package, so we refer to the [PlantGraphs.jl](https://github.com/VirtualPlantLab/PlantGraphs.jl) docs for specifics.
"""

# ╔═╡ 5f36f5e6-0245-4d61-aa19-239329af4fae
md"#### Rules"

# ╔═╡ 54ad8f3b-0a4f-449e-8f03-69e6dbab0400
md"""
The outermost stem segment grows another thinner stem segment and a branch every `branching_frequency` steps, if the stem is old enough to form a crown:
"""

# ╔═╡ e9c997bb-3daa-4546-b83e-0e7c534359ab
branching_rule = Rule( 
    StemTip,
    lhs = st -> data(st).age >= crown_start_age && data(st).age % 
		branching_frequency == 0,
    rhs = st -> Stem(data(st).D, data(st).height) + (
        BranchTip([branch_radius(data(st).D[1], Δr), segment_length], 
				  data(st).height, 0, needle_areas[data(st).branch_nr]),
        StemTip(data(st).D - [Δr, 0.0], data(st).height + segment_length, 
				data(st).branch_nr + 1, data(st).age + 1)
    )
);

# ╔═╡ 842a9074-3991-4284-bbf6-68cc37577bd8
md"The outermost stem segment grows another stem segment if not branching:"

# ╔═╡ a3575ab2-7ac7-4d5b-b705-9d2eee7eec0e
vertical_growth_rule = Rule( 
    StemTip,
    lhs = st -> !(data(st).age >= crown_start_age && data(st).age % 
		branching_frequency == 0),
    rhs = st -> Stem(data(st).D, data(st).height) + 
		StemTip(data(st).D, data(st).height + segment_length, 
				data(st).branch_nr, data(st).age + 1)
);

# ╔═╡ 6c9f3f1c-1d41-4540-8b27-83a457928d2e
md"The outermost branch segments grown another branch segment every `branchgrowth_frequency` steps:"

# ╔═╡ ac3335b1-088c-4663-82b2-53f4fe0f1893
branchgrowth_rule = Rule(
    BranchTip,
    lhs = br -> data(br).age % branchgrowth_frequency == 0,
    rhs = br -> Branch(data(br).D, data(br).height) +
        BranchTip(data(br).D, data(br).height, data(br).age + 1, 
				  data(br).needle_area)
);

# ╔═╡ 3a397c17-bc2e-4b30-a48b-a63585d984ef
md"The outermost branch segments also age if they are not growing:"

# ╔═╡ f93b1b41-5b37-485c-bfc0-92e0390bd60a
branchage_rule = Rule(
    BranchTip,
    lhs = br -> !(data(br).age % branchgrowth_frequency == 0),
    rhs = br -> BranchTip(data(br).D, data(br).height, data(br).age + 1, 
						  data(br).needle_area)
);

# ╔═╡ ccb00a00-398a-4b8f-bb7e-f35e5abf86e7
md"#### Instantiation & rewriting"

# ╔═╡ e4380613-b5ee-4eb8-9aab-7b15df42925a
md"Finally, we define our initial graph and apply the rewriting rules. For the final step, we remove the `StemTip` node, as we don't require it."

# ╔═╡ 00a7dc82-7a5a-4953-9d6f-4bbbdca181bd
axiom = StemTip([sapwood_r0, segment_length], segment_length/2, 1, 0);

# ╔═╡ 999107ec-0d08-484f-a9cf-5073590eb65f
begin
	plant = Graph(axiom = axiom, rules = 
	    (vertical_growth_rule, branching_rule, branchage_rule, branchgrowth_rule)
	)
	for _ in 1:(n_segments)
	    rewrite!(plant)
	end
	PlantGraphs.prune!(plant.graph, getid([node for node in getnodes(plant) if getstructmod(node) == :StemTip][1])) # remove StemTip
end

# ╔═╡ c383e15e-2adb-4dfa-af86-6980da48aa23
plotstructure(plant)

# ╔═╡ 0b74c50b-6fc8-4d50-9424-581968108fc7
md"### Connecting to the environment"

# ╔═╡ a5cbaceb-3387-4431-b4bd-b26867406721
md"""
We will again be using a single soil and air compartment, so there is little new in connecting the plant graph to the environment. However, for the visualisation of the resulting structure, removing the names of the structural modules with `names = ""` can give a clearer visualisation.
"""

# ╔═╡ 58001eea-6986-495d-8278-88080bf09a18
struct Soil <: PlantGraphs.Node end

# ╔═╡ 2cd061a1-0625-456a-9055-7de04f614384
struct Air <: PlantGraphs.Node end

# ╔═╡ 5830341b-92c7-45e9-9dcc-aca44570e667
graphs = [plant, Soil(), Air()];

# ╔═╡ b1bed12d-a8db-452c-9863-920344829a7e
intergraph_connections = [(1, 2) => (getnodes(plant)[1], :Soil), (1, 3) => (:BranchTip, :Air)];

# ╔═╡ 47a61b32-92e0-4100-a4fd-7a37c08d38a8
plantstructure = PlantStructure(graphs, intergraph_connections);

# ╔═╡ ad16819d-1b9e-45c3-8aa7-cca15f63ced7
plotstructure(plantstructure, names = "")

# ╔═╡ 6741bdad-518a-45e9-b127-863c6e9a36c5
md"## Function definition"

# ╔═╡ 58859042-be76-4e9f-a564-5d75d75d466b
md"### New functional modules"

# ╔═╡ 275b2d3b-9f36-4897-8c2a-328b4c2781e2
md"""
As mentioned before, the model has a set transpiration from the branch tips to the air, defined per unit needle area, and a total needle area for every branch tip.
To implement this, we need two new modules:
- A module for the branch tips describing their total needle area.
- A connection module between branch tips and the air that describes the transpiration in function of a given time-series for transpiration (per unit needle area) and the needle area of the branch tip segment.

Functionality is defined in terms of [`ModelingToolkit.jl`](https://docs.sciml.ai/ModelingToolkit/stable/) Systems, so we recommend unfamiliar users to look at their docs, specifically the [non-DSL way of defining systems](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/programmatically_generating/#The-Non-DSL-(non-@mtkmodel)-Way-of-Defining-an-System).
"""

# ╔═╡ b592211a-4db5-4797-869a-e358960a9836
md"We need our independent variable and derivative operator to be consistent with those used by the other functional modules, so we import then from the package first:"

# ╔═╡ 521ab69d-1346-427e-8f3f-66f869c7c568
import PlantModules: t, d

# ╔═╡ 8f280922-e866-4d2e-8ec1-6f59fa7d5102
md"#### Needle area module"

# ╔═╡ 6853391b-7018-4d47-a7f3-0199d65c56dd
md"The needle area of a branch is considered to be constant through time, so the functional module is pretty straightforward. The only tricky part is that we need to define it as a constant-valued variable rather than a parameter. This is because we need to reference this needle area in the transpiration connection, and we can only reference variables, not parameters."

# ╔═╡ 5bde3460-f174-4277-af4d-e46025cb9298
function needle_area_module(; name, needle_area)
    @parameters begin
        needle_area_val = needle_area, 
			[description = "Value of total needle area on branch", unit = u"cm^2"]
    end 
    @variables begin
        needle_area(t), [description = "Total needle area on branch", unit = u"cm^2"]
    end
    eqs = [needle_area ~ needle_area_val]
    return System(eqs, t; name)
end;

# ╔═╡ 6478d8cb-68ff-441a-981e-17d3f49875b6
md"#### Hydraulics module with effects of gravity"

# ╔═╡ 0a7ac994-eb69-400b-b904-73a7dbbc4b8a
function gravitational_hydraulic_module(; name, shape, ϕ_D, ϵ_D, Γ, T, D, P_0, M, height)
    D, ϕ_D, ϵ_D = 
		[PlantModules.correctdimensionality(shape, var) for var in [D, ϕ_D, ϵ_D]] 
        # turns scalar values into vectors of correct length

    num_D = getdimensionality(shape)
    R = 8.314
    ρ_w = 1.0
    g = 9.8 * 1e-5
    P = P_0 + R*T*M - ρ_w * g * height

    @constants (
        R = R, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"],
        P_unit = 1.0, [description = "Dummy constant for correcting units",
					   unit = u"MPa"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / cm^3"],
        g = g, [description = "Gravitational acceleration", unit = u"hN / g"] 
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϕ_D[1:num_D] = ϕ_D, [description = "Dimensional extensibility", 
							 unit = u"MPa^-1 * hr^-1"],
        ϵ_D[1:num_D] = ϵ_D, [description = "Dimensional elastic modulus", 
							 unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
        height = height, [description = "Height above ground level", unit = u"cm"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        Pₕ(t), [description = "Gravitational water potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", 
			   unit = u"mol / cm^3"],
        W(t) = volume(shape, D) / ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", 
				unit = u"MPa / hr", guess = 0.0],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment",
						 unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P + Π + Pₕ,
        Π ~ -R*T*M,
        Pₕ ~ ρ_w * g * height,
        ΔW ~ ΣF,
        V ~ W / ρ_w,
        V ~ volume(shape, D),
        [ΔD[i] ~ D[i]*ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100) +
			D[i]*ΔP/ϵ_D[i] for i in eachindex(D)]...,

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return System(eqs, t; name)
end;

# ╔═╡ bd152224-7f9b-4d1f-8649-550552fbf032
md"#### Fixed transpiration"

# ╔═╡ c18de4ad-dfcf-49de-9905-9b158ae0eae4
md"Next is the transpiration connection. We have time-series data for the transpiration per unit needle area, which we read in first. To use this in a functional module, we need to turn it into a function that will return the transpiration at _any_ timepoint, so some interpolation is required. We use a simple LinearInterpolation from the [`DataInterpolations.jl`](https://docs.sciml.ai/DataInterpolations/stable/) package for this."

# ╔═╡ 3b765f93-ed77-467b-a5ba-91d98f4a233a
readcsv(filepath) = readlines(filepath) .|>
	x -> split(x, ", ") .|> x -> parse(Float64, x);

# ╔═╡ deca4956-03d8-417a-8752-0a87fd6e26a0
transpiration_data = readcsv("./data/transpiration_data.csv");

# ╔═╡ 730ec309-f5c9-487f-b7b7-5fe47d8fb496
transpiration_rate = LinearInterpolation(last.(transpiration_data), first.(transpiration_data));

# ╔═╡ b7e20dd2-3f7f-42b3-af34-f4ee12cfe37b
plot(transpiration_rate, xlabel = "time (hr)", ylabel = "Specific transpiration (mg / s / m^2")

# ╔═╡ d0479f27-4c54-4561-83b5-a43ad1557f01
md"""
The actual module definition of the fixed transpiration connection is rather complicated. We have two main complications: Firstly, this specific connection has no inherent directionality, so we need to define it manually. Secondly, we need to relate the variables used in the connection to the variables of the nodes it connects.
"""

# ╔═╡ 914cba27-068f-4fa1-bc5d-5d7dbaa2645a
md"""
First, the specification of directionality. Let's compare with the usual hydraulic connection, for which water flow is defined in function of the water potential difference of the nodes it connects. Take the example of a connection between leaves and the air. If we consider the connection from the leaves' perspective or the air's perspective, we will both times have the same water potential difference but with a different sign. From both perspectives we will therefore also have the same water flow but with a different sign. This is perfect because our water flow from the air to the leaves needs to be the exact opposite of the flow from the leaves to the air. Because the flow direction is determined by the water potentials of the nodes it connects, we don't have to specify it ourselves.

Our fixed transpiration rate function, in contrast, does not give any indication what direction the water flow needs to go in, only the value. That is why we need to manually specify this when we define the connection. We can do this by using the special `original_order` keyword argument, which has the value `true` when the order between the structural modules of the node and the neighboring node is the same as originally specified by the user in the `connecting_modules` dictionary, and `false` otherwise. We'll here choose the order "branch first, air second". In this case, the water flow needs to be negative when the nodes are in the original order (water flows away from the tree during transpiration), and positive when not.
"""

# ╔═╡ 8ee56176-a3a0-447f-a275-36e337a325a5
md"""
The second complication is inherent to all connection modules: we need to specify how the variables from our connection (i.e. the graph edge) relate to those of the structural modules it connects (i.e. the graph nodes). This is done by defining a function with the following inputs:
- `node_MTK`: The functional module of the first node in the connection.
- `nb_node_MTK`: The functional module of the second node in the connection.
- `connection_MTK`: The functional module of the connection itself.
- `original_order` (optional): Whether the first and second nodes are in the original order specified in the `connecting_modules` dictionary.

This function is then returned as the second output of our connection function. For our model, the variable that needs to be connected is the needle area. As it is only defined for the `BranchTip` nodes and not the `Air` node, we need to use the `original_order` argument to make sure the variable is connected to the system of the `BranchTip` node.

For completeness' sake: The water flux variable `F` also needs to be connected to the node systems. However, the way we usually want to connect water fluxes to the nodes is by defining the net water influx of a node as the sum of *all* connected water fluxes. Because the variable of a node is connected to the variables of _multiple_ edge systems, this specification needs to be given separately. We don't need to change how water fluxes relate to the node systems for this model, but it is possible. See the `PlantCoupling` docstring for more details.
"""

# ╔═╡ 0c9d64be-93b1-4305-8de3-2291e065c303
function fixed_transpiration_connection(; name, original_order)
    @constants begin
        t_unit = 1, [description = "Dummy constant for correcting units", 
					 unit = u"hr"]
        uc = 1e-3 * 3600 * (1e-2)^2, [description = "Unit conversion from
			(mg / s / m^2) to (g / hr / cm^2)", unit = u"g / hr / cm^2"]
    end
    @variables begin
        F_s(t), [description = "Specific water flux from compartment 2 
			to compartment 1", unit = u"g / hr / cm^2"]
        F(t), [description = "Water flux from compartment 2 to compartment 1", 
			   unit = u"g / hr"]
        needle_area(t), [description = "Total needle area on branch", unit = u"cm^2"]
    end

    polarity = original_order ? -1 : 1

    eqs = [
        F_s ~ uc * transpiration_rate(t / t_unit),
        F ~ polarity * F_s * needle_area,
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, original_order) = (
        original_order ? 
			[connection_MTK.needle_area ~ node_MTK.needle_area] : [connection_MTK.needle_area ~ nb_node_MTK.needle_area]
    )

    return System(eqs, t; name), get_connection_eqset
end;

# ╔═╡ 3750e2ee-bec0-4f06-ba0f-4efb1e8b05f9
md"### Coupling"

# ╔═╡ 33f5bf30-fbac-4b96-bd3d-6aeb100540dd
md"""
The coupling of function and structure remains, luckily, straightforward. We only have to pay special attention to the order of the `BranchTip` and `Air` in `connecting_modules`, as this will determine how the `original_order` variable works for our connection module.
"""

# ╔═╡ 9abaf023-2ce5-4246-a41e-f3d507c77194
module_coupling = Dict(
	:Stem => [gravitational_hydraulic_module, constant_carbon_module, K_module],
    :Branch => [gravitational_hydraulic_module, constant_carbon_module, K_module],
	:BranchTip => [gravitational_hydraulic_module, needle_area_module,
				   constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module],
);

# ╔═╡ d81913d7-657d-4b81-a83c-e7185d6164bc
connecting_modules = Dict(
	(:Soil, :Stem) => constant_hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Branch) => hydraulic_connection,
    (:Branch, :Branch) => hydraulic_connection,
	(:Branch, :BranchTip) => hydraulic_connection,
	(:BranchTip, :Air) => fixed_transpiration_connection,
);

# ╔═╡ a411768f-fb22-4575-9666-b42cc0473d4c
plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ╔═╡ d26b79a0-5331-4714-8110-2e2574f1ef35
md"### Parameters"

# ╔═╡ e4a0ca52-3679-4d76-b63a-2d9aeaf31e41
md"""
The parameter specification mainly comes down to changing the default parameter values to those prescribed by the paper. Note that when we define new parameters or variables they need to be given a module-wide default value, even if those values are never used (here the case for `needle_area` and `height`). Some striking parameters are the zeros for the extensibility `ϕ_D` and metabolite concentration `M`: this was chosen because the model only considers elastic diameter changes and no irreversible growth, and does not consider carbon dynamics.
"""

# ╔═╡ 2ed83b9f-820a-4c0e-9ce4-42fcb35daade
begin
    default_changes = Dict(
		:needle_area => 0.0, :height => 0.0,
		:ϵ_D => ϵ_D_stem, :K_s => K_s_stem,
        :ϕ_D => 0.0, :M => 0.0, 
		:P_0 => PlantModules.soilfunc(0.5), 
		:shape => HollowCylinder(radial_fraction_sapwood)
	)

	module_defaults = Dict(
        :Soil => Dict(:W_max => 1e6, :W_r => 0.5),
    )

    connection_values = Dict(
        (:Soil, :Stem) => Dict(:K => K_roots),
    )

    plantparams = PlantParameters(; default_changes, module_defaults,
								  connection_values)
end;

# ╔═╡ c7f8bc7f-2fcd-4400-a681-699dce7541df
md"## Running the system"

# ╔═╡ 5025d163-dfc3-4a3f-b1eb-3d902f1d2c98
md"""
System generation and running also remains largely the same, with the exception that we specify to `solve` what solver we want. `FBDF` is a good solver for large, stiff systems of DAEs, as can be found on [`DifferentialEquations.jl`'s overview of ODE solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
"""

# ╔═╡ 4ee50611-6032-476d-ad2f-9393ad0b2200
system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false);

# ╔═╡ 7f99ec81-41c9-45bf-9e7d-2ab7dada677d
prob = ODEProblem(system, [], (0.0, 24.0), sparse = true);

# ╔═╡ b24016d3-cd1d-4ba8-b7a6-2ee8b41db1f1
sol = solve(prob, FBDF());

# ╔═╡ 087fa5db-83a0-4a47-b52a-7ddc6259a281
md"## Results"

# ╔═╡ 81aff1dd-ba69-498e-aa9c-775b76b42d2f
md"""
We now validate our model by comparing it to the simulation results and measured outputs found in the paper. Firstly, we verify that our transpiration was implemented correctly by plotting the total water influx of the air. After converting it to the units used in the paper, we can visually compare the results. Secondly, we compare our simulated water pressures for three different stem segments at different heights to those found in the paper. Finally, we compare our simulated diameter variation through time with the measured diameter variation from the paper.
"""

# ╔═╡ 305e53bb-f09f-4546-8865-1d2f6f74cae0
md"### Check transpiration"

# ╔═╡ 0368f643-ed99-475c-ac81-c17f9a82c8d3
md"""
When we either want to plot very specific variables, or if we want to use `ModelingToolkit.jl`'s plotting functionality to apply a transformation to variables before visualizing them, the `get_subsystem_variables` function can be used to extract the `Symbolic` representation of specific variables.
"""

# ╔═╡ eefbf375-b604-4a67-8861-0a428b088ee3
air_water_inflow = get_subsystem_variables(system, plantstructure, :ΣF, :Air)[1]

# ╔═╡ 38841aaf-d87b-48f2-9b5a-d800828a8791
md"""
We can then define a function that changes the transpiration from the units used in our model (g / hr) to those in the paper (mg / m^2 / s).
"""

# ╔═╡ 426b4aac-a610-488c-820f-198172c1e8b9
transpiration_unit_conversion(var) = var * 1e3 / (total_needle_area * (1e-2)^2) / 3600

# ╔═╡ 9e8bd19a-2790-41dd-8655-0120271c1d39
plot(sol, idxs = [transpiration_unit_conversion(air_water_inflow)], label = false,
	xlabel = "Time of day (h)", ylabel = "Transpiration (mg / m^2 / s)",
	xticks = 0:3:24, ylims = (0, 10), yticks = 0:2:10, color = :black)

# ╔═╡ e2ce9fbb-ff2a-4894-98e3-308d8b93c392
md"### Water tension"

# ╔═╡ 0accbe57-ab1c-435c-9812-0a173decb850
md"""
We use the same approach to select the pressure variable from `Stem` nodes at specific heights of the tree.
"""

# ╔═╡ b3186620-808b-481e-89f6-a037223ed990
begin
	stem_pressures = get_subsystem_variables(system, plantstructure, :P, :Stem)
	pressure_segment_heights = [10.0, crown_base_height, 1200.0]
	pressure_segment_nrs = ceil.(Int64, pressure_segment_heights / segment_length)
	
	plot(sol, idxs = stem_pressures[pressure_segment_nrs], ylims = (-0.4, 0.0),
		 xticks = 0.0:3.0:24.0, xlabel = "Time of day (h)", color = :black,
		 linewidth = [1 2 3], ylabel = "Water tension (MPa)",
		 label = ["Base of stem" "Crown base" "Top of tree"])
end

# ╔═╡ 28184995-7163-40ee-a773-eb414ecf2f28
md"### Diameter variation"

# ╔═╡ 59c6e1f7-ce34-4fcb-aa7b-7a650dc0a4d1
md"""
Finally, we can compare our simulation with actual measurements.
"""

# ╔═╡ 6b2ad036-b6f1-4753-a3f0-4a6d4dc6b7fd
begin
    dimension_variables = get_subsystem_variables(system, plantstructure, :D, :Stem)

    diameter_change_mm(var, sol) = 2*10*(var - sol[var][1]) 
		# substract first value and change from radius in cm to diameter in mm
	diameter_segment_heights = [6.5 * 1e2, 2.5 * 1e2]
    diameter_segment_nrs = ceil.(Int64, diameter_segment_heights ./ segment_length)
	
    diameter_datas = [readcsv("./data/diameter_data_$(idx).csv") 
					  for idx in ("c", "d")]
	
	subplots = [plot(), plot()]
	for i in eachindex(subplots)
		plot!(subplots[i], sol, idxs = 
			  [diameter_change_mm(
				  dimension_variables[diameter_segment_nrs[i]][1], sol
			  )],
			  label = "Simulated", color = :black, linewidth = 1)
		plot!(subplots[i], first.(diameter_datas[i]), last.(diameter_datas[i]), 
			  label = "Data", color = :black, linewidth = 2)
		plot!(subplots[i], 
			  title = "Stem at height $(diameter_segment_heights[i]) cm")
	end
	plot(subplots..., layout = (2, 1), ylabel = "Diameter change (mm)", 
		 size = (800, 600), margins = 5*Plots.mm, ylims = (-0.07, 0.01), 
		 yticks = -0.07:0.01:0.01, xticks = 0:3:24, xlabel = "Time of day (h)")
	savefig("fig_plantmodules_ex2_results.pdf")
end

# ╔═╡ 70df3d1e-82ad-4365-97ab-0e31283058c2
md"## Uncertainty analysis"

# ╔═╡ 9771ad64-1821-4f3e-bf63-6176ad7dce79
ϵ_D_r_range = [0.03 * 1e3, 0.27 * 1e3]

# ╔═╡ 97356ed2-d3cc-4007-b2d2-9e326be6c90c
sample_range(a, b) = (b-a) * rand() + a

# ╔═╡ d99adee7-c13a-44b2-8d8b-7e3d396751e3
begin
	uncertainty_plot = plot(ylabel = "Diameter change (mm)", 
		 size = (800, 600), margins = 5*Plots.mm, ylims = (-0.12, 0.01), 
		 yticks = -0.12:0.01:0.01, xticks = 0:3:24, xlabel = "Time of day (h)")
	for _ in 1:20
	    ϵ_D_r_sample = sample_range(ϵ_D_r_range...)
	    ϵ_D_sample = [ϵ_D_r_sample, 17.5 * ϵ_D_r_sample] #! 17.5 to variable
	
	    prob2 = remake_graphsystem(prob, system, plantstructure, :ϵ_D, 
								   [:Stem, :Branch, :BranchTip], ϵ_D_sample)
	    sol2 = solve(prob2, FBDF(), reltol = 1e-1)
	    plot!(uncertainty_plot, sol2, idxs = 		
			  [diameter_change_mm(dimension_variables[diameter_segment_nrs[2]][1], sol)], label = false, line_z = ϵ_D_r_sample)
	end
	uncertainty_plot
end

# ╔═╡ Cell order:
# ╟─e04d4d44-3795-49d0-91d3-645ff3c8265e
# ╟─61d4d897-e1ad-4f04-b0be-1df409c4e69f
# ╟─21271cd8-eef1-485a-978d-f46afdb753b6
# ╠═3fc5da40-acde-4398-8a53-99794e42ba94
# ╠═04361839-2e2f-4475-8830-277c4b07e3f6
# ╠═45991fe2-1dae-4633-b3c2-81b3e076c1eb
# ╠═8adf74c7-7fe0-42fd-bde5-9942f30fea36
# ╠═df79e96c-a021-4a24-ae8c-4090d014c1f8
# ╟─b2f4ec18-5447-4999-9da6-52ca770a8120
# ╟─ac068794-a6cc-4563-899a-4b5ca1b2c22b
# ╟─78f2ca6c-9a4e-4b3f-b0f5-d7df2a78fa36
# ╟─5a0bf882-68e2-4ad3-ae1f-43cac5d8f890
# ╠═1371c204-ba0a-4f18-9df3-d2614277bc11
# ╟─fdedf540-9c04-476f-9ffb-d7bc82e81f79
# ╟─e26238dd-d323-4161-9e95-bece1af1b3c5
# ╟─97b3b0b1-0408-46ba-96e0-01eadab25ea0
# ╟─479acd0d-ce44-411d-97bf-e12a49004a5d
# ╠═2a4b2b06-a057-4fef-b8fe-30700b39b564
# ╠═6d1a8b68-e8d3-43d9-84a3-8b3fdd600913
# ╠═1ca717b5-7c67-4727-a984-7fd6914a3257
# ╠═c0ce423e-ddec-4a2f-b724-8d55d3d05073
# ╠═17a090c0-f673-4c2d-b1e4-8f9f8a360179
# ╟─d59a64b2-9fdf-4723-b764-b27330a050fc
# ╟─c7150b44-966b-4c61-b66d-b8080cd51acd
# ╠═d7ea529f-4fe2-4daf-9c18-bb4621f0c86d
# ╠═ee63840c-392c-4e4a-88d9-a41ccb88fa6f
# ╠═dc48821b-ffb3-4599-b744-b69364e7ad6f
# ╠═578a128c-37e7-4b4d-8ba1-513565027110
# ╟─e319323b-6d26-410b-8574-411686d4067c
# ╟─b7d2469c-642e-435d-a399-5406c55e01c4
# ╟─5f36f5e6-0245-4d61-aa19-239329af4fae
# ╟─54ad8f3b-0a4f-449e-8f03-69e6dbab0400
# ╠═e9c997bb-3daa-4546-b83e-0e7c534359ab
# ╟─842a9074-3991-4284-bbf6-68cc37577bd8
# ╠═a3575ab2-7ac7-4d5b-b705-9d2eee7eec0e
# ╟─6c9f3f1c-1d41-4540-8b27-83a457928d2e
# ╠═ac3335b1-088c-4663-82b2-53f4fe0f1893
# ╟─3a397c17-bc2e-4b30-a48b-a63585d984ef
# ╠═f93b1b41-5b37-485c-bfc0-92e0390bd60a
# ╟─ccb00a00-398a-4b8f-bb7e-f35e5abf86e7
# ╟─e4380613-b5ee-4eb8-9aab-7b15df42925a
# ╠═00a7dc82-7a5a-4953-9d6f-4bbbdca181bd
# ╠═999107ec-0d08-484f-a9cf-5073590eb65f
# ╠═c383e15e-2adb-4dfa-af86-6980da48aa23
# ╟─0b74c50b-6fc8-4d50-9424-581968108fc7
# ╟─a5cbaceb-3387-4431-b4bd-b26867406721
# ╠═58001eea-6986-495d-8278-88080bf09a18
# ╠═2cd061a1-0625-456a-9055-7de04f614384
# ╠═5830341b-92c7-45e9-9dcc-aca44570e667
# ╠═b1bed12d-a8db-452c-9863-920344829a7e
# ╠═47a61b32-92e0-4100-a4fd-7a37c08d38a8
# ╠═ad16819d-1b9e-45c3-8aa7-cca15f63ced7
# ╟─6741bdad-518a-45e9-b127-863c6e9a36c5
# ╟─58859042-be76-4e9f-a564-5d75d75d466b
# ╟─275b2d3b-9f36-4897-8c2a-328b4c2781e2
# ╟─b592211a-4db5-4797-869a-e358960a9836
# ╠═521ab69d-1346-427e-8f3f-66f869c7c568
# ╟─8f280922-e866-4d2e-8ec1-6f59fa7d5102
# ╟─6853391b-7018-4d47-a7f3-0199d65c56dd
# ╠═5bde3460-f174-4277-af4d-e46025cb9298
# ╟─6478d8cb-68ff-441a-981e-17d3f49875b6
# ╠═0a7ac994-eb69-400b-b904-73a7dbbc4b8a
# ╟─bd152224-7f9b-4d1f-8649-550552fbf032
# ╟─c18de4ad-dfcf-49de-9905-9b158ae0eae4
# ╠═3b765f93-ed77-467b-a5ba-91d98f4a233a
# ╠═deca4956-03d8-417a-8752-0a87fd6e26a0
# ╠═730ec309-f5c9-487f-b7b7-5fe47d8fb496
# ╠═b7e20dd2-3f7f-42b3-af34-f4ee12cfe37b
# ╟─d0479f27-4c54-4561-83b5-a43ad1557f01
# ╟─914cba27-068f-4fa1-bc5d-5d7dbaa2645a
# ╟─8ee56176-a3a0-447f-a275-36e337a325a5
# ╠═0c9d64be-93b1-4305-8de3-2291e065c303
# ╟─3750e2ee-bec0-4f06-ba0f-4efb1e8b05f9
# ╟─33f5bf30-fbac-4b96-bd3d-6aeb100540dd
# ╠═9abaf023-2ce5-4246-a41e-f3d507c77194
# ╠═d81913d7-657d-4b81-a83c-e7185d6164bc
# ╠═a411768f-fb22-4575-9666-b42cc0473d4c
# ╟─d26b79a0-5331-4714-8110-2e2574f1ef35
# ╟─e4a0ca52-3679-4d76-b63a-2d9aeaf31e41
# ╠═2ed83b9f-820a-4c0e-9ce4-42fcb35daade
# ╟─c7f8bc7f-2fcd-4400-a681-699dce7541df
# ╟─5025d163-dfc3-4a3f-b1eb-3d902f1d2c98
# ╠═4ee50611-6032-476d-ad2f-9393ad0b2200
# ╠═7f99ec81-41c9-45bf-9e7d-2ab7dada677d
# ╠═b24016d3-cd1d-4ba8-b7a6-2ee8b41db1f1
# ╟─087fa5db-83a0-4a47-b52a-7ddc6259a281
# ╟─81aff1dd-ba69-498e-aa9c-775b76b42d2f
# ╟─305e53bb-f09f-4546-8865-1d2f6f74cae0
# ╟─0368f643-ed99-475c-ac81-c17f9a82c8d3
# ╠═eefbf375-b604-4a67-8861-0a428b088ee3
# ╟─38841aaf-d87b-48f2-9b5a-d800828a8791
# ╠═426b4aac-a610-488c-820f-198172c1e8b9
# ╠═9e8bd19a-2790-41dd-8655-0120271c1d39
# ╟─e2ce9fbb-ff2a-4894-98e3-308d8b93c392
# ╟─0accbe57-ab1c-435c-9812-0a173decb850
# ╠═b3186620-808b-481e-89f6-a037223ed990
# ╟─28184995-7163-40ee-a773-eb414ecf2f28
# ╟─59c6e1f7-ce34-4fcb-aa7b-7a650dc0a4d1
# ╠═6b2ad036-b6f1-4753-a3f0-4a6d4dc6b7fd
# ╟─70df3d1e-82ad-4365-97ab-0e31283058c2
# ╠═9771ad64-1821-4f3e-bf63-6176ad7dce79
# ╠═97356ed2-d3cc-4007-b2d2-9e326be6c90c
# ╠═d99adee7-c13a-44b2-8d8b-7e3d396751e3
