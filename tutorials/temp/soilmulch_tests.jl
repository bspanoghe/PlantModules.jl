using Revise, Infiltrator #!
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
import PlantGraphs: Node
using ModelingToolkit, OrdinaryDiffEq, Unitful, Plots
using SoilMulch

# # Structure
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

graphs = [plant_graph, Soil(), Air()]

intergraph_connections = [[1, 2] => (:Root, :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)] #! order matters now

struct_connections = PlantStructure(graphs, intergraph_connections)

# # Function

# ## SoilMulch

num_days = 10
pr = (sh = 0.15, sw = 0.25, sstar = 0.55, ks = 160, β = 16.0, n = 0.50, zr = 30, cws = 0.20, ce = 1.0,
    ct = 1.0, μ = 0.11, ϕ = 0.60, zm = 5.00, ϑh = 0.05, km = 230, g = 6.50, tsen = 90, γ = 0.005, cmax = 0.98,
    r = 0.39, dt = 1/24
);
rain = rainfall_poisson(num_days/pr.dt, 0.5, 0.9*pr.dt);
swb = sol_swb_crop(rain, 0.5, 0.1, 0.35, pr)
plot(scatter(rain, ylabel = "Rain (cm?)", label = false), plot(swb.s, ylabel = "Soil moisture (vol fraction?)", label = false))




C_root = 300e-6 # (mol/cm^3 = 10^6 µmol/g) We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 400e-6 # Idem here
C_leaf = 450e-6 # And here!

default_changes = Dict(:Γ => 0.4, :P => 0.2, :T => 293.15)

module_defaults = Dict(
	:Root => Dict(:shape => Cuboid([5.0, 10.0, 20.0], [0.7, 0.1, 0.05]), :D => [30, 3, 1], :M => C_root), #! changed C to M_const #! include check whether module_defaults variabkle names correspond with default_params / default_u0s ?
	:Stem => Dict(:shape => Cylinder([6.0, 15.0], [0.8, 0.03]), :D => [1.5, 10], :M => C_stem),
	:Leaf => Dict(:shape => Sphere([3.0], [0.45]), :M => C_leaf),
	:Soil => Dict(:W_max => 666.0, :T => 288.15, :K => 100.0),
	:Air => Dict(:W_max => 1e5, :W_r => 0.7, :K => 0.1)
)

connecting_modules = [
	(:Soil, :Root) => (hydraulic_connection, Dict()),
	(:Root, :Stem) => (hydraulic_connection, Dict()), # 6*10^-7 kg/s/MPa * 1000 g/kg * 3600 s/h = 2.2 g/h/MPa
	(:Stem, :Leaf) => (hydraulic_connection, Dict()),
	(:Leaf, :Air) => (hydraulic_connection, Dict()),
	(:Soil, :Air) => (const_hydraulic_connection, Dict(:K => 0.01)) #! check value
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
prob = ODEProblem(sys_simpl, [], (0.0, 5*24), sparse = true)
@time sol = solve(prob);

plotgraph(sol, graphs[1]) .|> display;
plotgraph(sol, graphs[2]) .|> display;
plotgraph(sol, graphs[3]) .|> display;
plotgraph(sol, graphs[1:3], varname = :Ψ)
plotgraph(sol, graphs[1:3], varname = :W)

plotgraph(sol, graphs[1], varname = [:P])
plotnode(sol, PlantModules.getnodes(graphs[1])[1], varname = [:P, :V]) .|> display