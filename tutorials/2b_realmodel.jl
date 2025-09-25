using Revise
using Plots
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful

# # Structure definition

# ## Preparatory calculations

# ### stem segments and change in radius
segment_length = 20 # input
tree_length = 1200 # input
stem_r0 = 12.0 # input
branching_frequency = 3 # dealer's choice

n_segments = tree_length ÷ segment_length
Δr = stem_r0/n_segments

# ### crown
crown_length = tree_length - 600 # input
n_crown_segments = crown_length ÷ segment_length
segments_before_crown = (n_segments - n_crown_segments)

# ### branch radius
function branch_radius(r, Δr)
    area_0 = r^2*pi
    area_new = (r - Δr)^2*pi
    area_branch = area_0 - area_new  # from da vinci's rule (total cross area remains constant)
    return sqrt(area_branch / pi)
end

# ### total needle length

get_cylinder_length(area, radius) = (area - 2*radius^2*pi) / (2*radius*pi)
    # A = 2 * (D[1]^2 * pi) + (2 * D[1] * pi) * D[2]
    # D[2] = (A - 2 * (D[1]^2 * pi)) /  (2 * D[1] * pi)

needle_area = 51 * 100^2 # input
needle_radius = 0.06 # based on https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.13465
total_needle_length = get_cylinder_length(needle_area, needle_radius)

# ### sanity check: estimate number of needles in tree assuming 10 cm needles
num_needles = needle_area / (needle_radius^2 * pi * 10)

# ### needles
n_needle_segments = n_crown_segments ÷ 3
needle_lengths = [i * total_needle_length / sum(1:n_needle_segments) for i in n_needle_segments:-1:1]

# ### root dimensions

rootarea = 96 * 100^2 # input
rootradius = 0.015 # based on https://cdnsciencepub.com/doi/abs/10.1139/x98-206
total_root_length = get_cylinder_length(rootarea, rootradius) #! hm.


# ## Actual definition

Base.@kwdef mutable struct StemTip <: Node
    tsb::Integer = 1 # time since branching
    nd::Integer = 1 # number of divisions / branchings
    D::Vector = [stem_r0, segment_length]
end

Base.@kwdef mutable struct Stem <: Node
    D::Vector
end

Base.@kwdef mutable struct Branch <: Node
    D::Vector
    age::Integer = 1
end

Base.@kwdef mutable struct Needles <: Node
    D::Vector
end

Base.@kwdef mutable struct Roots <: Node
    D::Vector = [rootradius, total_root_length] 
end

# ## Growing rules

# ### Plant parts

stop_rule = Rule(
    StemTip,
    lhs = gt -> data(gt).D[1] < Δr
)

vertical_growth_rule = Rule(
    StemTip,
    lhs = gt -> data(gt).tsb < (segments_before_crown+branching_frequency) && data(gt).D[1] > Δr,
    rhs = gt -> Stem(data(gt).D) + StemTip(data(gt).tsb + 1, data(gt).nd, data(gt).D - [Δr, 0.0])
)

branching_rule = Rule(
    StemTip,
    lhs = gt -> data(gt).tsb == (segments_before_crown+branching_frequency),
    rhs = gt -> Stem(data(gt).D) + (
        Branch([branch_radius(data(gt).D[1], Δr), segment_length], 0) + Needles([needle_radius, needle_lengths[data(gt).nd]]),
        StemTip(segments_before_crown+1, data(gt).nd + 1, data(gt).D - [Δr, 0.0])
    )
)

branchgrowth_rule = Rule(
    Branch,
    lhs = br -> (data(br).age == 1 && length(children(br)) == 1),
    rhs = br -> Branch(data(br).D, data(br).age + 1) + Branch(data(br).D, 0)
)

# ## Instantiation

axiom = Roots() + StemTip()
plant = Graph(axiom = axiom, rules = (stop_rule, vertical_growth_rule, branching_rule, branchgrowth_rule))
for _ in 1:(n_segments+10)
    rewrite!(plant)
end
PlantStructure([plant], []) |> plotstructure

[node for node in getnodes(plant) if PlantModules.getstructmod(node) == :Stem]