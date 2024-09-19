using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs
# using MultiScaleTreeGraph
using ModelingToolkit, DifferentialEquations, Unitful
using Plots; import GLMakie.draw
using BenchmarkTools

# Structure

## Plant

mutable struct Stem <: Node
	D::Vector
end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

function branchingrule(x)
    leaf_D = [2.0, 1.0, 0.1]
    parent_D = data(x).D
    child_D = [parent_D[1] / sqrt(2), 1.0]
    return Stem(parent_D) + (Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)), Stem(child_D) + (Leaf(leaf_D), Leaf(leaf_D)))
end

function branchlessrule(x)
    parent_D = data(x).D
    return Stem(parent_D) + Stem(parent_D)
end

function growify!(plant, growthrate)
    for node in apply(plant, Query(Stem))
        node.D = [sqrt(growthrate)*node.D[1], growthrate*node.D[2]]
    end

	for node in apply(plant, Query(Leaf))
        node.D = [growthrate*node.D[1], growthrate*node.D[2], node.D[3]]
    end
end

rule1 = Rule(Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchingrule
)

rule2 = Rule(Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchlessrule
)

axiom = Stem([0.5, 5.0]) + (Leaf([3.0, 1.0, 0.1]), Leaf([5.0, 3.0, 0.1]))

plant = Graph(axiom = axiom, rules = (rule1,))
num_iterations = 5
for _ in 1:num_iterations
    rewrite!(plant)
	growify!(plant, 1.1)
end

# draw(plant)
convert_to_MTG(plant) |> PlantModules.MultiScaleTreeGraph.DataFrame

## Environment

soil_graph = Soil()
air_graph = Air()

## Combined

graphs = [plant, soil_graph, air_graph]
intergraph_connections = [[1, 2] => (PlantModules.root(plant), :Soil), [1, 3] => (:Leaf, :Air), [2, 3] => (:Soil, :Air)]
struct_connections = [graphs, intergraph_connections]


# Functional processes

@independent_variables t, [description = "Time", unit = u"hr"]; #! add documentation
d = Differential(t);

function simpler_hydraulic_module(; name, T, shape::PlantModules.Shape, Γ, P, D)
    num_D = length(shape.ϵ_D)
    @constants (
        P_0 = 0.0, [description = "Minimum pressure", unit = u"MPa"],
        R = 8.314, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϵ_D[1:num_D] = shape.ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        ϕ_D[1:num_D] = shape.ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        Γ = Γ, [description = "Critical turgor pressure", unit = u"MPa"],
        ρ_w = 1.0, [description = "Density of water", unit = u"g / cm^3"],
        MPa_unit = 1.0, [description = "Dummy variable for correcting units", unit = u"MPa"], #!
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation () #! extend validation function so L is ok?
        W(t) = PlantModules.volume(shape, D) * ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t) = PlantModules.volume(shape, D), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P - Π, # Water potential consists of a solute- and a pressure component
        Π ~ R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        V ~ W / ρ_w, # Volume is directly related to water content  
        V ~ PlantModules.volume(shape, D), # Volume is also directly related to compartment dimensions
        [ΔD[i] ~ D[i] * ϕ_D[i] * max(P - Γ, P_0) for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure

        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    # return ODESystem(eqs, t; name, continuous_events = [P ~ Γ]) #!
    return ODESystem(eqs, t; name)
end


default_params = merge(PlantModules.default_params, 
	(simpler_hydraulic_module = (T = 298.15, shape = Sphere(ϵ_D = [1.0], ϕ_D = [1.0]), Γ = 0.3),),
)
default_u0s = merge(PlantModules.default_u0s,
	(simpler_hydraulic_module = (P = 0.5, D = [15],),),
)

C_root = 300e-6
C_stem = 400e-6
C_leaf = 450e-6

module_defaults = (
	Stem = (shape = PlantModules.Cilinder(ϵ_D = [0.6, 1.5], ϕ_D = 1e-5 .* [8, 3]), D = [1.5, 10], M = C_stem),
	Leaf = (shape = PlantModules.Cuboid(ϵ_D = [0.3, 0.3, 10.0], ϕ_D = 1e-5 .* [3, 3, 0.1]), M = C_leaf),
	Soil = (W_max = 10000.0, T = 288.15),
	Air = (W_r = 0.8,)
)

connecting_modules = [
	(:Soil, :Stem) => (PlantModules.hydraulic_connection, [:K => 80]),
	(:Stem, :Stem) => (PlantModules.hydraulic_connection, [:K => 80]),
	(:Stem, :Leaf) => (PlantModules.hydraulic_connection, [:K => 40]),
	(:Leaf, :Air) => (PlantModules.hydraulic_connection, [:K => 1e-3]),
	(:Soil, :Air) => (PlantModules.hydraulic_connection, [:K => 5e-3]) #! check value
] # values based on https://www.mdpi.com/2073-4441/10/8/1036

func_connections = [connecting_modules, PlantModules.multi_connection_eqs]

# Coupling 

module_coupling = [
	PlantModules.hydraulic_module => [:Stem, :Leaf],
	PlantModules.constant_carbon_module => [:Stem, :Leaf],
	PlantModules.environmental_module => [:Soil, :Air],
	PlantModules.Ψ_soil_module => [:Soil],
	PlantModules.Ψ_air_module => [:Air]
]

# Rev her up

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24))
@time sol = solve(prob);

# Plotting
plotgraph(sol, graphs[1], func_varname = :W)
plotgraph(sol, graphs[1], func_varname = :W, ylims = (0, 10))

plotgraph(sol, graphs[1], func_varname = :P)
plotgraph(sol, graphs[1], func_varname = :M)

plotgraph(sol, graphs[2], func_varname = :W)
plotgraph(sol, graphs[1], func_varname = :D)
plotgraph(sol, graphs[1:2], func_varname = :Ψ)
plotgraph(sol, graphs[1], func_varname = :ΣF)

plotnode(sol, PlantModules.root(plant), func_varname = :W)
plotnode(sol, PlantModules.nodes(plant)[end], func_varname = :D)






function sizedep_hydraulic_connection(; name, K, connection_shape)
    @parameters (
        K_s = K_s, [description = "Specific hydraulic conductivity of connection", unit = u"g / hr / MPa / cm^2"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
		D_1(t), [description = "Dimensions of compartment 1", unit = u"cm^2"],
		D_2(t), [description = "Dimensions of compartment 2", unit = u"cm^2"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1)
		K ~ K_s * (PlantModules.surface_area(connection_shape, D_1) + PlantModules.surface_area(connection_shape, D_2))/2
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order) = [ 
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
		connection_MTK.D_1 ~ node_MTK.D,
        connection_MTK.D_2 ~ nb_node_MTK.D,
    ]

    return ODESystem(eqs, t; name), get_connection_eqset
end







import ModelingToolkit: get_unknowns, get_name

get_MTKunknown_symbol(s::SymbolicUtils.Symbolic) = s.metadata[ModelingToolkit.VariableSource][2]

## Pry the ODE system corresponding with given node out of the ODE solution
function getnodesystem(sol::ODESolution, node)
    system = getfield(sol.prob.f, :sys)
    nodename = string(PlantModules.structmod(node)) * string(PlantModules.id(node))
    nodesystem = [subsys for subsys in getproperty(system, :systems) if get_name(subsys) == Symbol(nodename)][1]
    return nodesystem
end

nodesystem = getnodesystem(sol, PlantModules.root(graphs[1]))
func_varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystem)] |> unique

@which sol[getproperty(nodesystem, :Ψ)]

for func_varname in func_varnames
    println(func_varname); @time sol[getproperty(nodesystem, func_varname)]
end