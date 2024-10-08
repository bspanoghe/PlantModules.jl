# # Tutorial 2: A real model

# ## Introduction

# In this tutorial, we'll cover more advanced functionality of the package with the goal of making a more practically usable plant model.
using Revise #!
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, MultiScaleTreeGraph
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots; import GLMakie.draw

# ## Structure

# ### The plant
# For this model, we'll assume we have a **file** containing the plant's structure, which we may have gotten from some other plant modeling program.
# Currently, the package only has a function for reading XEG files. However, any file format can be used provided the user is able to convert it to a graph.

plant_graph = readXEG("./tutorials/temp/structures/beech10.xeg") #! change to beech

# The graph still requires some processing to make it suitable for modeling with PlantModules.
# Luckily, the [MultiScaleTreeGraph.jl](https://github.com/VEZY/MultiScaleTreeGraph.jl) package has some excellent functionality for processing graphs.
# We can use the `convert_to_MTG` function to change our graph into this package's format allowing us to use all its functions.
mtg = convert_to_MTG(plant_graph)

# Currently the nodes have their dimensions defined in separate variables. The hydraulic module we intend to use expects a vector `D` instead, so we need to change this.
DataFrame(mtg, [:diameter, :length, :width])

function combine_dimensions(l, d, w)
	if all(isnothing.([l, d, w]))
		return nothing
	elseif isnothing(w) # no width defined => shape is a cilinder
		return 1e2*[d/2, l] # 1e2 to go from m to cm
	else # otherwise we're dealing with a cuboid
		return 1e2*[l, w, 5e-4] # leaf thickness assumed to be 0.5 mm
	end
end

transform!(mtg, [:length, :diameter, :width] => combine_dimensions => :D)
DataFrame(mtg, [:D]) # inspect the results

# Rename the "ShortShoot" nodes to "Shoot" for brevity's sake
traverse!(mtg, node -> symbol!(node, "Shoot"), symbol = "ShortShoot")
# Delete nodes without any dimensions defined. These are for graph visualisation purposes but are not actual physical structures.
mtg = delete_nodes!(mtg, filter_fun = node -> isnothing(node.D))

if length(mtg) > 100 #!
	toomuch = [node for node in PlantModules.nodes(mtg[1][1][1])]
	mtg = delete_nodes!(mtg, filter_fun = node -> node in toomuch)
end
DataFrame(mtg)

convert_to_PG(mtg) |> draw


# ### Environment
# This part is the same as previous tutorial.

struct Soil <: PlantGraphs.Node end
struct Air <: PlantGraphs.Node end

soil_graph = Soil()
air_graph = Air()

graphs = [mtg, soil_graph, air_graph]
intergraph_connections = [[1, 2] => (mtg, :Soil), [1, 3] => (:Leaf, :Air)]
struct_connections = PlantStructure(graphs, intergraph_connections)

# ## Function

# ### Defining new functional modules

# Previous tutorial, we simply made use of the pre-implemented functional modules. Now, let's take a look at how we can define our own. 
# We'll start by making a simple model for carbon dynamics.

# First thing we need is to define how much photosynthetically active radiation (PAR) our leaves receive. 
# We'd like to define this as a general Julia function rather than a differential equation so that we have more freedom in our definition.
# This allows you to, for example, write a function which interpolates between observed PAR values. 
# We'll just use a sine function bound to the positive values to simulate a simple day-night cycle:

get_PAR_flux(t) = max(0, 400 * sin(t/24*2*pi - 8))

# We also need the actual photosynthesis model.
# We could define this ourselves using differential equations, or simply use the existing implementation from [PlantBiophysics.jl](https://github.com/VEZY/PlantBiophysics.jl).
# All we have to do is wrap this model in a julia function with the desired inputs and outputs and we're done.

using PlantBiophysics, PlantBiophysics.PlantMeteo, PlantSimEngine
using Memoization

# PlantBiophysics models are fast, but using memoization will still speed up the final model up significantly.
# Memoization essentially makes the function remember what the outputs are for inputs it has seen before, rather than having to compute them again.

@memoize function get_assimilation_rate(PAR_flux, T, LAI, k)
	Kelvin_to_C = -273.15
	meteo = Atmosphere(T = T + Kelvin_to_C, Wind = 1.0, P = 101.3, Rh = 0.65, Ri_PAR_f = PAR_flux)
	m = ModelList(
		Fvcb(), # calculate CO2 assimilation rate
		Medlyn(0.03, 0.92), # calculate stomatal conductance, see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2486.2010.02375.x
		Beer(k), # calculate amount of light intercepted
		status = (Tₗ = meteo[:T], LAI = LAI, Cₛ = meteo[:Cₐ], Dₗ = meteo[:VPD], RI_PAR_f = meteo[:Ri_PAR_f])
	)
	run!(m, meteo)
	return only(m[:A]) |> x -> max(x, 0) # extract result of the first (and only) timestep
end

# Most (complex) functions need to be registered using `@register_symbolic` before we can use them in ModelingToolkit.
@register_symbolic get_assimilation_rate(PAR_flux, T, LAI, k)

import PlantModules: t, d #!

# Finally we define the actual photosynthesis module. The one defined here is very simple for illustration purposes.
function photosynthesis_module(; name, T, M, shape)
	@constants (
		uc1 = (10^-6 * 10^-4 * 60^2), [description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)", unit = u"(mol/cm^2/hr) / (µmol/m^2/s)"],
		# the output from PlantBiophysics.jl is in different units than we use for our ODEs, so we need to change this
	)
	@parameters (
		T = T, [description = "Temperature", unit = u"K"],
		LAI = 2.0, [description = "Leaf Area Index", unit = u"cm^2 / cm^2"],
		k = 0.5, [description = "Light extinction coefficient", unit = u"N/N"],
		carbon_decay_rate = 0.3, [description = "Rate at which carbon is consumed for growth", unit = u"hr^-1"],
	)
	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
		PF(t), [description = "Incoming PAR flux", unit = u"J / s / m^2"], #! make sure not to use variable name from other func mod used in same struct mod
		A(t), [description = "Carbon assimilation rate", unit = u"µmol / m^2 / s"],
		D(t)[1:length(shape.ϵ_D)], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		PF ~ get_PAR_flux(t)
		A ~ get_assimilation_rate(PF, T, LAI, k)
        d(M) ~ uc1 * A * cross_area(shape, D) / volume(shape, D) - carbon_decay_rate*M
    ]
    return ODESystem(eqs, t; name, checks = false) #! checks back to true?
end

# This part is the same as in the first tutorial, except that we'll also make use of the non-constant `hydraulic_connection`.
# For this connection, the hydraulic conductance is related to the cross-area of the plant part.
# The way we defined the leaf shape, its hydraulic conductance will be proportional to its surface area.
# This is why this hydraulic conductance should only be used for a connection with the air, and connections with plant parts need to use a constant hydraulic connection.
module_defaults = Dict(
	:Internode => Dict(:shape => Cilinder(ϵ_D = [2.0, 8.0], ϕ_D = 1e-3 * [2, 2]), :M => 300e-6),
	:Shoot => Dict(:shape => Cilinder(ϵ_D = [2.0, 8.0], ϕ_D = 1e-3 * [3, 3]), :M => 350e-6),
	:Leaf => Dict(:shape => Cuboid(ϵ_D = [1.5, 1.5, 10.0], ϕ_D = 1e-3 * [3, 3, 0.1]), :M => 400e-6, :K_s => 5e-4),
	:Soil => Dict(:W_max => 1e9, :T => 293.15), #! W_max
	:Air => Dict(:K => 1e-1)
)

connecting_modules = [
	(:Soil, :Internode) => (const_hydraulic_connection, Dict()),
    (:Internode, :Internode) => (hydraulic_connection, Dict()),
	(:Internode, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Shoot) => (hydraulic_connection, Dict()),
	(:Shoot, :Leaf) => (const_hydraulic_connection, Dict()),
	(:Internode, :Leaf) => (const_hydraulic_connection, Dict()),
    (:Leaf, :Air) => (hydraulic_connection, Dict())
]

func_connections = PlantFunctionality(module_defaults = module_defaults, connecting_modules = connecting_modules)

# ## Coupling

module_coupling = Dict(
	:Internode => [hydraulic_module, constant_carbon_module, sizedep_K_module],
	:Shoot => [hydraulic_module, constant_carbon_module, sizedep_K_module],
	:Leaf => [hydraulic_module, constant_carbon_module, sizedep_K_module],
	:Soil => [environmental_module, Ψ_soil_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
)

# ## Generate and run system

system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false)

sys_simpl = structural_simplify(system)
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
@time sol = solve(prob);

# ## Plotting

plotgraph(sol, graphs[1], varname = :W, structmod = :Leaf)
plotgraph(sol, graphs[1], varname = :M, structmod = :Leaf)
plotgraph(sol, graphs[1], varname = :W, ylims = (0, 0.1))
plotgraph(sol, graphs[2], varname = :W)

plotgraph(sol, graphs[1:2], varname = :Ψ)