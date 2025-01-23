using Revise, Infiltrator #!
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs, MultiScaleTreeGraph
import PlantGraphs: Node
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots
using BenchmarkTools

mutable struct Root <: Node end
mutable struct Stem <: Node end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

if !isdefined(Main, :plant_graph) # for debugging mostly
	plant_graph = Root() + Stem() + (Leaf([2]), Leaf([2.5]))
end

# draw(plant_graph, resolution = (500, 400))

soil_graph = Soil()
air_graph = Air()

graphs = [plant_graph, soil_graph, air_graph]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = PlantStructure(graphs, intergraph_connections)


C_root = 300e-6 # (mol/cm^3 = 10^6 µmol/g) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 400e-6 # Idem here
C_leaf = 450e-6 # And here!
# loosely based on https://www.researchgate.net/figure/Total-soluble-sugar-content-in-the-leaf-a-and-root-tissues-b-of-Homjan-HJ_fig3_226164026

PlantModules.default_values

default_changes = Dict(:Γ => 0.4, :P => 0.2, :T => 293.15)

module_defaults = Dict(
	:Root => Dict(:shape => Cuboid([5.0, 10.0, 20.0], [0.7, 0.1, 0.05]), :D => [30, 3, 1], :M => C_root), #! changed C to M_const #! include check whether module_defaults variabkle names correspond with default_params / default_u0s ?
	:Stem => Dict(:shape => Cylinder([6.0, 15.0], [0.8, 0.03]), :D => [1.5, 10], :M => C_stem),
	:Leaf => Dict(:shape => Sphere([3.0], [0.45]), :M => C_leaf),
	:Soil => Dict(:W_max => 500.0, :T => 288.15),
	:Air => Dict(:W_r => 0.8, :K => 0.001)
)

connecting_modules = [
	(:Soil, :Root) => (hydraulic_connection, Dict()),
	(:Root, :Stem) => (hydraulic_connection, Dict()), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Stem, :Leaf) => (hydraulic_connection, Dict()),
	(:Leaf, :Air) => (hydraulic_connection, Dict()),
	(:Soil, :Air) => (hydraulic_connection, Dict()) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

func_connections = PlantFunctionality(; default_changes, module_defaults, connecting_modules)

module_coupling = Dict(
	:Root => [hydraulic_module, constant_carbon_module, K_module],
	:Stem => [hydraulic_module, constant_carbon_module, K_module],
	:Leaf => [hydraulic_module, constant_carbon_module, K_module],
	:Soil => [environmental_module, Ψ_soil_module, constant_K_module],
	:Air => [environmental_module, Ψ_air_module, constant_K_module],
)

system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, [], (0.0, 5*24))
@time sol = solve(prob);

plotgraph(sol, graphs[1]) .|> display
plotgraph(sol, graphs[2]) .|> display
plotgraph(sol, graphs[3]) .|> display
plotgraph(sol, graphs[1:2], varname = :Ψ)
plotgraph(sol, graphs[1:3], varname = :W)

plotgraph(sol, graphs[1], varname = [:D, :V], structmod = [:Stem, :Leaf]) .|> display
plotnode(sol, PlantModules.getnodes(graphs[1])[1], varname = [:P, :V]) .|> display