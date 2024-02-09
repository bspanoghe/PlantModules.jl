using Pkg; Pkg.activate(".")
include("PlantModules.jl")
using PlantGraphs, ModelingToolkit, GLMakie, DifferentialEquations #! MTK imports etc. should not be necessary when package is done

mutable struct Root <: Node end
mutable struct Stem <: Node end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end


C_root = 0.5 # We'll assume the soluble carbon content remains constant over the simulated time
C_stem = 1 # Idem here
C_leaf = 3 # And here!
Ψ_soil_func(W_r) = -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19) # An empirical relationship between the soil water potential and relative water content
Ψ_air_func(T, W_r) = R * T / V_w * log(W_r) #! What's this equation called again? (Ask Jeroen)

@constants R = 8.314 V_w = 18e-6

PlantModules.default_vals

model_defaults = [:Γ => 0.4, :P => 0.2, :M => 15.0, :T => 293.15]
module_defaults = [
	:Root => [:D => [0.3, 0.05, 0.03], :ϵ_D => [5.0, 0.3, 0.2], :ϕ_D => [0.7, 0.1, 0.05], :C => C_root],
	:Stem => [:D => [0.4, 0.03], :ϵ_D => [6.0, 0.15], :ϕ_D => [0.8, 0.03], :C => C_stem],
	:Leaf => [:ϵ_D => [3.0], :ϕ_D => [0.45], :C => C_leaf],
	:Soil => [:W_max => 500.0, :T => 288.15, :Ψ => Ψ_soil_func],
	:Air => [:Ψ => Ψ_air_func]
]

module_coupling = [
	PlantModules.hydraulic_module => [:Root, :Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air]
]

plant_graph = Root() + Stem() + (Leaf([0.25]), Leaf([0.15]))
# draw(plant_graph, resolution = (500, 400))

soil_graph = Soil()
air_graph = Air()

intergraph_connections = [(:Air, :Leaf), (:Soil, :Root)]

struct_connections = [[plant_graph, soil_graph, air_graph], intergraph_connections]

func_connections = [
	:default => PlantModules.hydraulic_connection,
	(:Soil, :Root) => [:K => 50],
	(:Root, :Stem) => [:K => 800],
	(:Leaf, :Air) => [:K => 1e-3]
]

plantsys = PlantModules.PlantSystem(
	model_defaults = model_defaults,
	module_defaults = module_defaults,
	module_coupling = module_coupling,
	struct_connections = struct_connections,
	func_connections = func_connections
)

time_span = (0, 7*24.0) # We'll simulate our problem for a timespan of 7 days
prob = ODEProblem(plantsys, time_span)
sol = solve(prob)

plot(sol, struct_modules = [:Soil], func_vars = [:W]) #! imagine that this works



#######

graphs = struct_connections[1]
for graph in graphs
	MTK_systems = []
    MTK_connections = []
    MTK_u0 = []
	iteratedescendants(graph, 
		add_MTK_info!(node, 
			model_defaults = model_defaults, module_defaults = module_defaults,
			module_coupling = module_coupling, struct_connections = struct_connections,
			func_connections = func_connections, MTK_systems = MTK_systems,
			MTK_connections = MTK_connections, MTK_u0 = MTK_u0))
end

# Apply a function to a node and all its descendants
function iteratedescendants(node, graph, func::Function; kwargs...)
	func(node; kwargs...)
	for chnode in PlantModules.children(node, graph)
		iteratedescendants(chnode, graph, func; kwargs...)
	end
end

# Default behaviour: start from graph root
iteratedescendants(graph, func::Function; kwargs...) = iteratedescendants(PlantModules.root(graph), graph, func; kwargs...)

function add_MTK_info!(node; model_defaults, module_defaults, module_coupling, struct_connections,
	func_connections, MTK_systems, MTK_connections, MTK_u0)

	getMTKsystem(node)

end

function getMTKsystem(node, module_coupling)

	structmodule = PlantModules.nodetype(node)
	func_module = [coupling.first for coupling in module_coupling if structmodule in coupling.second][1]
	MTKsystem = func_module(; name = string(structmodule))

end

phew(x; a = 3, b = 4) = x+a+b
dump(:(phew(x, a = 5, b = 100)))


# tests #

## iteratedescendants 
mutable struct Foo <: Node
	bar::Int
end
testgraph = Foo(1) + (Foo(2), Foo(3) + (Foo(7), Foo(9)))
bars = Int[]
iteratedescendants(testgraph, (x; extra) -> push!(bars, x.data.bar + extra), extra = 3)
bars == [1, 2, 3, 7, 9] .+ 3