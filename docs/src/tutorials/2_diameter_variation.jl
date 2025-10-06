### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ 3fc5da40-acde-4398-8a53-99794e42ba94
using Pkg; Pkg.activate("../..")

# ╔═╡ 7081348a-c47b-4e0c-ae4c-cafe8ccb2b4e
using Plots

# ╔═╡ 04361839-2e2f-4475-8830-277c4b07e3f6
using PlantModules

# ╔═╡ 45991fe2-1dae-4633-b3c2-81b3e076c1eb
using PlantGraphs

# ╔═╡ 9915cee8-eb00-4b48-9f72-80f6ef7e68d0
using ModelingToolkit, OrdinaryDiffEq, Unitful

# ╔═╡ 8adf74c7-7fe0-42fd-bde5-9942f30fea36
using DataInterpolations

# ╔═╡ 78f2ca6c-9a4e-4b3f-b0f5-d7df2a78fa36
md"## Preparatory calculations"

# ╔═╡ 1371c204-ba0a-4f18-9df3-d2614277bc11
begin
    # ## dimensions of compartments

    # ### stem segments
    segment_length = 0.2 * 1e2 # cm; from m (input)
    tree_length = 13.2 * 1e2 # cm; from m (input)
    n_segments = tree_length ÷ segment_length

    # ### stem radius at tree base
    stem_r0 = 0.2 * 1e2 / 2 # cm; from the diameter in m (input)
    radial_fraction_sapwood = 0.30 # (estimate)
    sapwood_r0 = stem_r0 * radial_fraction_sapwood

    # ### crown
    crown_base_height = 6.0 * 1e2 # cm; from m (input)
    crown_length = tree_length - crown_base_height
    n_crown_segments = crown_length ÷ segment_length
    crown_start_age = (n_segments - n_crown_segments)

    # ### branching
    branching_frequency = 3 # (estimate)
    n_branches = n_crown_segments ÷ branching_frequency

    Δr = sapwood_r0/(n_branches + 1) # (+1 to prevent stem radius from hitting 0)

    # ### branch radius
    # calculate using da vinci's rule: total cross area remains constant before and after branching
    function branch_radius(r, Δr)
        area_0 = r^2*pi
        area_new = (r - Δr)^2*pi
        area_branch = area_0 - area_new
        return sqrt(area_branch / pi)
    end

    # ### needle areas
    total_needle_area = 50.8 * 1e2^2 # cm^2; from m^2 (input)
    needle_areas = [i * total_needle_area / sum(1:n_branches) for i in n_branches:-1:1]

    # ### branchgrowth frequency
    branchgrowth_frequency = 3 # how many steps it takes branches to grow a new segment (estimate)

    # ### root dimensions
    total_root_area = 96.1 * 1e2^2 # cm^2; from m^2 (input)

    # ## parameters

    # ### hydraulic conductivities
    permeability = 2.5e-12 * 1e2^2 # cm^2; from m^2 (input)
    η = 1.0 * 1e-9 * (1/3600) # MPa hr; from mPa s (dynamic viscosity of water)
    l = segment_length # cm
    ρ_w = 1.0 # g / cm^3 (density of water)
    K_s_stem = ρ_w * permeability / η / l # g / hr / MPa / cm^2

    L_p = 3.2e-8 * 1e2 * 3600 # cm / hr / MPa; from m / s / MPa (input)
    K_roots = ρ_w * L_p * total_root_area # g / hr / MPa

    # ### elastic modulus
    ϵ_D_stem = [0.1 * 1e3, 17.5 * 0.1 * 1e3] # MPa; from GPa (input)
end

# ╔═╡ fdedf540-9c04-476f-9ffb-d7bc82e81f79
md"## Structure definition"

# ╔═╡ 97b3b0b1-0408-46ba-96e0-01eadab25ea0
md"### Defining a new shape"

# ╔═╡ 2a4b2b06-a057-4fef-b8fe-30700b39b564
struct HollowCylinder{T} <: PlantModules.Shape
    frac_sapwood::T # fraction of radius that is conducting sapwood
    function HollowCylinder(frac_sapwood::T) where {T <: Number}
        @assert 0 <= frac_sapwood <= 1
        return new{T}(frac_sapwood)
    end
end

# ╔═╡ 6d1a8b68-e8d3-43d9-84a3-8b3fdd600913
PlantModules.getdimensionality(::HollowCylinder) = 2

# ╔═╡ 1ca717b5-7c67-4727-a984-7fd6914a3257
PlantModules.cross_area(hc::HollowCylinder, D::AbstractArray) = D[1]^2 * pi * (2 / hc.frac_sapwood - 1)

# ╔═╡ c0ce423e-ddec-4a2f-b724-8d55d3d05073
PlantModules.volume(hc::HollowCylinder, D::AbstractArray) = cross_area(hc, D) * D[2]

# ╔═╡ d59a64b2-9fdf-4723-b764-b27330a050fc
md"## Structural module types"

# ╔═╡ dc48821b-ffb3-4599-b744-b69364e7ad6f
mutable struct StemTip <: Node
    D::Vector
    height
    branch_nr
    age
end

# ╔═╡ d7ea529f-4fe2-4daf-9c18-bb4621f0c86d
mutable struct Stem <: Node
    D::Vector
    height
end

# ╔═╡ 578a128c-37e7-4b4d-8ba1-513565027110
mutable struct BranchTip <: Node
    D::Vector
    height
    age
    needle_area
end

# ╔═╡ ee63840c-392c-4e4a-88d9-a41ccb88fa6f
mutable struct Branch <: Node
    D::Vector
    height
end

# ╔═╡ 5f36f5e6-0245-4d61-aa19-239329af4fae
md"### Growing rules"

# ╔═╡ e9c997bb-3daa-4546-b83e-0e7c534359ab
branching_rule = Rule( # top stem segment grows another thinner stem segment and a branch every `branching_frequency` steps, if stem is old enough to form crown
    StemTip,
    lhs = st -> data(st).age >= crown_start_age && data(st).age % branching_frequency == 0,
    rhs = st -> Stem(data(st).D, data(st).height) + (
        BranchTip([branch_radius(data(st).D[1], Δr), segment_length], data(st).height, 0, needle_areas[data(st).branch_nr]),
        StemTip(data(st).D - [Δr, 0.0], data(st).height + segment_length, data(st).branch_nr + 1, data(st).age + 1)
    )
)

# ╔═╡ a3575ab2-7ac7-4d5b-b705-9d2eee7eec0e
vertical_growth_rule = Rule( # top stem segment grows only another stem segment if not branching
    StemTip,
    lhs = st -> !(data(st).age >= crown_start_age && data(st).age % branching_frequency == 0),
    rhs = st -> Stem(data(st).D, data(st).height) + 
        StemTip(data(st).D, data(st).height + segment_length, data(st).branch_nr, data(st).age + 1)
)

# ╔═╡ ac3335b1-088c-4663-82b2-53f4fe0f1893
branchgrowth_rule = Rule( # top branch segments grown another branch segment every `branchgrowth_frequency` steps
    BranchTip,
    lhs = br -> data(br).age % branchgrowth_frequency == 0,
    rhs = br -> Branch(data(br).D, data(br).height) +
        BranchTip(data(br).D, data(br).height, data(br).age + 1, data(br).needle_area)
)

# ╔═╡ f93b1b41-5b37-485c-bfc0-92e0390bd60a
branchage_rule = Rule( # branch tips also age if not growing
    BranchTip,
    lhs = br -> !(data(br).age % branchgrowth_frequency == 0),
    rhs = br -> BranchTip(data(br).D, data(br).height, data(br).age + 1, data(br).needle_area)
)

# ╔═╡ e319323b-6d26-410b-8574-411686d4067c
md"## Instantiation"

# ╔═╡ 00a7dc82-7a5a-4953-9d6f-4bbbdca181bd
axiom = StemTip([sapwood_r0, segment_length], segment_length/2, 1, 0)

# ╔═╡ 999107ec-0d08-484f-a9cf-5073590eb65f
begin
	plant = Graph(axiom = axiom, rules = 
	    (vertical_growth_rule, branching_rule, branchage_rule, branchgrowth_rule)
	)
	for _ in 1:(n_segments)
	    rewrite!(plant)
	end
	PlantGraphs.prune!(plant.graph, getid([node for node in getnodes(plant) if getstructmod(node) == :StemTip][1]))
end

# ╔═╡ c383e15e-2adb-4dfa-af86-6980da48aa23
plotstructure(plant)

# ╔═╡ 0b74c50b-6fc8-4d50-9424-581968108fc7
md"### Environment and connections"

# ╔═╡ 58001eea-6986-495d-8278-88080bf09a18
struct Soil <: Node end

# ╔═╡ 2cd061a1-0625-456a-9055-7de04f614384
struct Air <: Node end

# ╔═╡ 5830341b-92c7-45e9-9dcc-aca44570e667
graphs = [plant, Soil(), Air()]

# ╔═╡ b1bed12d-a8db-452c-9863-920344829a7e
intergraph_connections = [(1, 2) => (getnodes(plant)[1], :Soil), (1, 3) => (:BranchTip, :Air)]

# ╔═╡ 47a61b32-92e0-4100-a4fd-7a37c08d38a8
plantstructure = PlantStructure(graphs, intergraph_connections)

# ╔═╡ ad16819d-1b9e-45c3-8aa7-cca15f63ced7
plotstructure(plantstructure, names = "")

# ╔═╡ 6741bdad-518a-45e9-b127-863c6e9a36c5
md"## Function"

# ╔═╡ 58859042-be76-4e9f-a564-5d75d75d466b
md"### New functional modules"

# ╔═╡ 521ab69d-1346-427e-8f3f-66f869c7c568
import PlantModules: t, d

# ╔═╡ 3b765f93-ed77-467b-a5ba-91d98f4a233a
readcsv(filepath) = readlines(filepath) .|>
	x -> split(x, ", ") .|> x -> parse(Float64, x)

# ╔═╡ 8f280922-e866-4d2e-8ec1-6f59fa7d5102
md"#### Needle area module"

# ╔═╡ 5bde3460-f174-4277-af4d-e46025cb9298
function needle_area_module(; name, needle_area)
    @parameters begin
        needle_area_val = needle_area, [description = "Value of total needle area on branch", unit = u"cm^2"]
    end 
    @variables begin
        needle_area(t), [description = "Total needle area on branch", unit = u"cm^2"]
    end
    eqs = [needle_area ~ needle_area_val]
    return System(eqs, t; name)
end

# ╔═╡ 6478d8cb-68ff-441a-981e-17d3f49875b6
md"#### Hydraulics module with effects of gravity"

# ╔═╡ 0a7ac994-eb69-400b-b904-73a7dbbc4b8a
function gravitational_hydraulic_module(; name, shape, ϕ_D, ϵ_D, Γ, T, D, P_0, M, height)
    D, ϕ_D, ϵ_D = [PlantModules.correctdimensionality(shape, var) for var in [D, ϕ_D, ϵ_D]] 
        # turns scalar values into vectors of correct length

    num_D = getdimensionality(shape)
    R = 8.314
    ρ_w = 1.0
    g = 9.8 * 1e-5
    P = P_0 + R*T*M - ρ_w * g * height

    @constants (
        R = R, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
        P_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / cm^3"],
        g = g, [description = "Gravitational acceleration", unit = u"hN / g"] # (from N / kg) Pa = N/m^2 => MPa = hN/cm^2
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϕ_D[1:num_D] = ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        ϵ_D[1:num_D] = ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
        height = height, [description = "Height above ground level", unit = u"cm"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        Pₕ(t), [description = "Gravitational water potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation ()
        W(t) = volume(shape, D) / ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr", guess = 0.0], #!
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P + Π + Pₕ,
        Π ~ -R*T*M,
        Pₕ ~ ρ_w * g * height,
        ΔW ~ ΣF,
        V ~ W / ρ_w,
        V ~ volume(shape, D),
        [ΔD[i] ~ D[i]*ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100) + D[i]*ΔP/ϵ_D[i] for i in eachindex(D)]...,

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return System(eqs, t; name)
end

# ### Fixed transpiration

# #### transpiration rate from data

# ╔═╡ deca4956-03d8-417a-8752-0a87fd6e26a0
transpiration_data = readcsv("./tutorials/clouddata/transpiration_data.csv")

# ╔═╡ 730ec309-f5c9-487f-b7b7-5fe47d8fb496
transpiration_rate = LinearInterpolation(last.(transpiration_data), first.(transpiration_data))

# ╔═╡ b7e20dd2-3f7f-42b3-af34-f4ee12cfe37b
plot(transpiration_rate) # in mg / s / m^2

# #### module

# ╔═╡ 0c9d64be-93b1-4305-8de3-2291e065c303
function fixed_transpiration_connection(; name)
    @constants begin
        t_unit = 1, [description = "Dummy constant for correcting units", unit = u"hr"]
        uc = 1e-3 * 3600 * (1e-2)^2, [description = "Unit conversion from (mg / s / m^2) to (g / hr / cm^2)", unit = u"g / hr / cm^2"]
    end
    @variables begin
        F_s(t), [description = "Specific water flux from compartment 2 to compartment 1", unit = u"g / hr / cm^2"]
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"]
        needle_area(t), [description = "Total needle area on branch", unit = u"cm^2"]
    end

    polarity = occursin("Air", split(string(name), "_")[1]) ? 1 : -1 #! ugly :( => use `correct_order` as input of function?

    eqs = [
        F_s ~ uc * transpiration_rate(t / t_unit),
        F ~ polarity * F_s * needle_area,
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, correct_order) = (
        correct_order ? [connection_MTK.needle_area ~ node_MTK.needle_area] : [connection_MTK.needle_area ~ nb_node_MTK.needle_area]
    )

    return System(eqs, t; name), get_connection_eqset
end


# ## Coupling

# ╔═╡ 9abaf023-2ce5-4246-a41e-f3d507c77194
module_coupling = Dict(
	:Stem => [gravitational_hydraulic_module, constant_carbon_module, K_module],
    :Branch => [gravitational_hydraulic_module, constant_carbon_module, K_module],
	:BranchTip => [gravitational_hydraulic_module, needle_area_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module],
)

# ╔═╡ d81913d7-657d-4b81-a83c-e7185d6164bc
connecting_modules = Dict(
	(:Soil, :Stem) => constant_hydraulic_connection,
	(:Stem, :Stem) => hydraulic_connection,
	(:Stem, :Branch) => hydraulic_connection,
    (:Branch, :Branch) => hydraulic_connection,
	(:Branch, :BranchTip) => hydraulic_connection,
	(:BranchTip, :Air) => fixed_transpiration_connection,
)

# ╔═╡ a411768f-fb22-4575-9666-b42cc0473d4c
plantcoupling = PlantCoupling(; module_coupling, connecting_modules)

# # Parameters

# ╔═╡ 2ed83b9f-820a-4c0e-9ce4-42fcb35daade
begin
    module_defaults = Dict(
        :Stem => Dict(:ϵ_D => ϵ_D_stem, :K_s => K_s_stem),
        :Branch => Dict(:ϵ_D => ϵ_D_stem, :K_s => K_s_stem),
        :BranchTip => Dict(:ϵ_D => ϵ_D_stem, :K_s => K_s_stem),
        :Soil => Dict(:W_max => 1e6, :W_r => 0.5),
    )

    default_changes = Dict(
        :needle_area => 0.0, :height => 0.0, :ϕ_D => 0.0, :M => 0.0, :P_0 => PlantModules.soilfunc(0.5),
        :shape => HollowCylinder(radial_fraction_sapwood))

    connection_values = Dict(
        (:Soil, :Stem) => Dict(:K => K_roots),
    )

    plantparams = PlantParameters(; default_changes, module_defaults, connection_values)
end

# # Run it

# ╔═╡ 4ee50611-6032-476d-ad2f-9393ad0b2200
system = generate_system(plantstructure, plantcoupling, plantparams, checkunits = false)

# ╔═╡ 7f99ec81-41c9-45bf-9e7d-2ab7dada677d
prob = ODEProblem(system, [], (0.0, 24.0), sparse = true)

# @time remake_graphsystem!(prob, system, plantstructure, :ϕ_D, :Stem, [10.0, 10.0])

# ╔═╡ b24016d3-cd1d-4ba8-b7a6-2ee8b41db1f1
@time sol = solve(prob, FBDF())

# # Show it

# ╔═╡ eb3f7b2e-805d-47b6-80ea-840de59e9905
plotgraph(sol, plantstructure, structmod = [:Soil, :Stem, :Branch], varname = :Ψ)

# ╔═╡ b3186620-808b-481e-89f6-a037223ed990
plotgraph(sol, plantstructure, structmod = [:Stem], varname = :P, xticks = 0.0:3.0:24.0)

# ╔═╡ 6b2ad036-b6f1-4753-a3f0-4a6d4dc6b7fd
begin
    dimension_variables = get_subsystem_variables(system, plantstructure, :D, :Stem)

    diameter_change_mm(var, sol) = 10*(var - sol[var][1])
    low_segment_nr = ceil(Int64, 2.5 / 0.2)

    low_diameter_data = readcsv("./tutorials/clouddata/diameter_data_d.csv")
    begin
        p_diameter = plot(ylabel = "Diameter change (mm)", size = (800, 400), margins = 5*Plots.mm)#, ylims = (-0.07, 0.01), yticks = -0.07:0.01:0.01)
        plot!(p_diameter, sol, idxs = [diameter_change_mm(dimension_variables[low_segment_nr][1], sol)], label = "Simulated")
        plot!(p_diameter, first.(low_diameter_data), last.(low_diameter_data), label = "Data")
    end
end

# ╔═╡ 10e4162a-8bd1-449e-8b85-f1f5670d8cee
plot(sol, idxs = [get_subsystem_variables(system, plantstructure, :F, (:Stem, :Stem))[2] / 1000],
    xlims = (6.0, 18.0), ylabel = "Sap flow (kg/h)", xlabel = "Time of day (h)", legend = false)

# ╔═╡ f193a527-f243-43e4-96e2-2046d619cc5b
plotgraph(sol, plantstructure, structmod = :Soil, varname = :ΣF)

# ╔═╡ b21f4bea-6190-436a-a160-8d761cea0653
plotgraph(sol, plantstructure, structmod = :Soil, varname = :Ψ)

# ╔═╡ 9e8bd19a-2790-41dd-8655-0120271c1d39
plotgraph(sol, plantstructure, structmod = :Air, varname = :ΣF)

# check transpiration

# ╔═╡ 6b13d3b3-cd71-4f2c-b96b-5b53c1e25cbd
transp_max_paper = 8.0 # mg / m^2 / s

# ╔═╡ ab2c22a9-05d0-4c17-9b79-dc5beadbf46c
transp_max_expected = transp_max_paper * 1e-3 * (total_needle_area * (1e-2)^2) * 3600 # g / hr

# ╔═╡ Cell order:
# ╠═7081348a-c47b-4e0c-ae4c-cafe8ccb2b4e
# ╠═3fc5da40-acde-4398-8a53-99794e42ba94
# ╠═04361839-2e2f-4475-8830-277c4b07e3f6
# ╠═45991fe2-1dae-4633-b3c2-81b3e076c1eb
# ╠═9915cee8-eb00-4b48-9f72-80f6ef7e68d0
# ╠═8adf74c7-7fe0-42fd-bde5-9942f30fea36
# ╟─78f2ca6c-9a4e-4b3f-b0f5-d7df2a78fa36
# ╠═1371c204-ba0a-4f18-9df3-d2614277bc11
# ╟─fdedf540-9c04-476f-9ffb-d7bc82e81f79
# ╟─97b3b0b1-0408-46ba-96e0-01eadab25ea0
# ╠═2a4b2b06-a057-4fef-b8fe-30700b39b564
# ╠═6d1a8b68-e8d3-43d9-84a3-8b3fdd600913
# ╠═1ca717b5-7c67-4727-a984-7fd6914a3257
# ╠═c0ce423e-ddec-4a2f-b724-8d55d3d05073
# ╟─d59a64b2-9fdf-4723-b764-b27330a050fc
# ╠═dc48821b-ffb3-4599-b744-b69364e7ad6f
# ╠═d7ea529f-4fe2-4daf-9c18-bb4621f0c86d
# ╠═578a128c-37e7-4b4d-8ba1-513565027110
# ╠═ee63840c-392c-4e4a-88d9-a41ccb88fa6f
# ╟─5f36f5e6-0245-4d61-aa19-239329af4fae
# ╠═e9c997bb-3daa-4546-b83e-0e7c534359ab
# ╠═a3575ab2-7ac7-4d5b-b705-9d2eee7eec0e
# ╠═ac3335b1-088c-4663-82b2-53f4fe0f1893
# ╠═f93b1b41-5b37-485c-bfc0-92e0390bd60a
# ╟─e319323b-6d26-410b-8574-411686d4067c
# ╠═00a7dc82-7a5a-4953-9d6f-4bbbdca181bd
# ╠═999107ec-0d08-484f-a9cf-5073590eb65f
# ╠═c383e15e-2adb-4dfa-af86-6980da48aa23
# ╟─0b74c50b-6fc8-4d50-9424-581968108fc7
# ╠═58001eea-6986-495d-8278-88080bf09a18
# ╠═2cd061a1-0625-456a-9055-7de04f614384
# ╠═5830341b-92c7-45e9-9dcc-aca44570e667
# ╠═b1bed12d-a8db-452c-9863-920344829a7e
# ╠═47a61b32-92e0-4100-a4fd-7a37c08d38a8
# ╠═ad16819d-1b9e-45c3-8aa7-cca15f63ced7
# ╟─6741bdad-518a-45e9-b127-863c6e9a36c5
# ╟─58859042-be76-4e9f-a564-5d75d75d466b
# ╠═521ab69d-1346-427e-8f3f-66f869c7c568
# ╠═3b765f93-ed77-467b-a5ba-91d98f4a233a
# ╟─8f280922-e866-4d2e-8ec1-6f59fa7d5102
# ╠═5bde3460-f174-4277-af4d-e46025cb9298
# ╠═6478d8cb-68ff-441a-981e-17d3f49875b6
# ╠═0a7ac994-eb69-400b-b904-73a7dbbc4b8a
# ╠═deca4956-03d8-417a-8752-0a87fd6e26a0
# ╠═730ec309-f5c9-487f-b7b7-5fe47d8fb496
# ╠═b7e20dd2-3f7f-42b3-af34-f4ee12cfe37b
# ╠═0c9d64be-93b1-4305-8de3-2291e065c303
# ╠═9abaf023-2ce5-4246-a41e-f3d507c77194
# ╠═d81913d7-657d-4b81-a83c-e7185d6164bc
# ╠═a411768f-fb22-4575-9666-b42cc0473d4c
# ╠═2ed83b9f-820a-4c0e-9ce4-42fcb35daade
# ╠═4ee50611-6032-476d-ad2f-9393ad0b2200
# ╠═7f99ec81-41c9-45bf-9e7d-2ab7dada677d
# ╠═b24016d3-cd1d-4ba8-b7a6-2ee8b41db1f1
# ╠═eb3f7b2e-805d-47b6-80ea-840de59e9905
# ╠═b3186620-808b-481e-89f6-a037223ed990
# ╠═6b2ad036-b6f1-4753-a3f0-4a6d4dc6b7fd
# ╠═10e4162a-8bd1-449e-8b85-f1f5670d8cee
# ╠═f193a527-f243-43e4-96e2-2046d619cc5b
# ╠═b21f4bea-6190-436a-a160-8d761cea0653
# ╠═9e8bd19a-2790-41dd-8655-0120271c1d39
# ╠═6b13d3b3-cd71-4f2c-b96b-5b53c1e25cbd
# ╠═ab2c22a9-05d0-4c17-9b79-dc5beadbf46c
