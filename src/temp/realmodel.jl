using Pkg; Pkg.activate(".")
include("../PlantModules.jl"); using .PlantModules
using PlantGraphs, ModelingToolkit, DifferentialEquations, Unitful, Plots, MultiScaleTreeGraph
import GLMakie.draw

# Creating plant #

plant_graph = readXEG("./src/temp/structures/test.xeg") #! change to beech
mtg = convert_to_MTG(plant_graph)

DataFrame(mtg, [:diameter, :length, :width])

function combine_dimensions(l, d, w)
	if all(isnothing.([l, d, w]))
		return nothing
	elseif isnothing(w)
		return [l, d]
	else
		return [l, w, 1.5e-4]
	end
end

transform!(mtg, [:length, :diameter, :width] => combine_dimensions => :D)

DataFrame(mtg, [:D])

me_structmods = [PlantModules.structmod(node) for node in PlantModules.nodes(mtg)] |> unique

for me_structmod in me_structmods
	num_nodes = length([node for node in PlantModules.nodes(mtg) if PlantModules.structmod(node) == me_structmod])
	println("There are $num_nodes nodes of type $me_structmod")
end

# E #

C_root = 250 # (mol/m^3) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 300 # Idem here
C_leaf = 330 # And here!
# loosely based on https://www.researchgate.net/figure/Total-soluble-sugar-content-in-the-leaf-a-and-root-tissues-b-of-Homjan-HJ_fig3_226164026

module_defaults = (
	Root = (shape = PlantModules.Cuboid(ϵ_D = [5.0, 0.3, 0.2], ϕ_D = [0.7, 0.1, 0.05]), D = [0.3, 0.03, 0.01], M = C_root),
	Pnode = (), #! it just aint right
    Internode = (shape = PlantModules.Cilinder(ϵ_D = [6.0, 0.15], ϕ_D = [0.8, 0.03]), D = [0.015, 0.1], M = C_stem),
    Meristem = (),
	Leaf = (shape = PlantModules.Sphere(ϵ_D = [3.0], ϕ_D = [0.45]), M = C_leaf),
	Soil = (W_max = 2000.0, T = 288.15),
	Air = ()
)

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]



soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = [graphs, intergraph_connections]

connecting_modules = [
	() => PlantModules.hydraulic_connection,
	(:Soil, :Root) => (PlantModules.hydraulic_connection, [:K => 8]),
	(:Root, :Stem) => (PlantModules.hydraulic_connection, [:K => 4]), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-2]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-2]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

get_connection_eqs = PlantModules.hydraulic_connection_eqs #!

func_connections = [connecting_modules, get_connection_eqs]

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))	 #! Generate warning for missing variable defaults that aren't dummies
sol = solve(prob)

PlantModules.plotgraph(sol, graphs[1]) .|> display
PlantModules.plotgraph(sol, graphs[2]) .|> display
PlantModules.plotgraph(sol, graphs[3]) .|> display