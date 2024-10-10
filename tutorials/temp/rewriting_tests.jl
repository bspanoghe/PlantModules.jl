using BenchmarkTools, Infiltrator, Revise
using Pkg; Pkg.activate("./tutorials")
using PlantModules
using PlantGraphs
using ModelingToolkit, OrdinaryDiffEq, Unitful
using Plots; import GLMakie.draw


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

rule = Rule(
    Stem,
    lhs = node -> !has_descendant(node, condition = n -> data(n) isa Stem)[1],
    rhs = branchingrule
)

axiom = Stem([0.5, 5.0]) + (Leaf([3.0, 1.0, 0.1]), Leaf([5.0, 3.0, 0.1]))

plant = Graph(axiom = axiom, rules = (rule,))
num_iterations = 3
for _ in 1:num_iterations
    rewrite!(plant)
end

convert_to_MTG(plant) |> PlantModules.MultiScaleTreeGraph.DataFrame

## Environment

soil_graph = Soil()
air_graph = Air()

## Combined

graphs = [plant, soil_graph, air_graph]
intergraph_connections = [(1, 2) => (PlantModules.root(plant), :Soil), (1, 3) => (:Leaf, :Air), (2, 3) => (:Soil, :Air)]
struct_connections = PlantStructure(graphs, intergraph_connections)

# Functional processes

module_defaults = Dict(
	:Stem => Dict(:shape => Cilinder(ϵ_D = [2.0, 4.5], ϕ_D = 1e-3 .* [8, 3]),
        :D => [1.5, 10], :M => 400e-6, :K_s => 10),
	:Leaf => Dict(:shape => Cuboid(ϵ_D = [1.5, 1.5, 10.0], ϕ_D = 1e-3 .* [3, 3, 0.1]),
        :M => 450e-6, :K_s => 1e-5),
	:Soil => Dict(:W_max => 10000.0, :T => 288.15, :K => 10),
	:Air => Dict(:W_r => 0.8, :K => 0.1)
)

connecting_modules = [
	(:Soil, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Stem) => (hydraulic_connection, Dict()),
	(:Stem, :Leaf) => (const_hydraulic_connection, Dict(:K => 10)),
	(:Leaf, :Air) => (hydraulic_connection, Dict()),
	(:Soil, :Air) => (hydraulic_connection, Dict())
]

func_connections = PlantFunctionality(module_defaults = module_defaults, connecting_modules = connecting_modules)

# Coupling 

module_coupling = Dict(
    :Stem => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Leaf => [hydraulic_module, constant_carbon_module, sizedep_K_module],
    :Soil => [environmental_module, Ψ_soil_module, constant_K_module],
    :Air => [environmental_module, Ψ_air_module, constant_K_module],
)

# Rev her up

system = generate_system(struct_connections, func_connections, module_coupling, checkunits = false);
sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 5*24));
@time sol = solve(prob);

wah = [sys.W for sys in system.systems if :W in get_MTKunknown_symbol.(sys.unknowns)]
@time wooh = sol[wah];

#=
SOLTIMES

3 (45): 0.1s
4: 1.8s
5 (189): 6s
6 (381): 28s
7 (765): 135s
=#

histogram(sol.t)

# Plotting

plotgraph(sol, graphs[1], varname = :W)
plotgraph(sol, graphs[2], varname = :W)

plotgraph(sol, graphs[1:2], varname = :Ψ)
plotgraph(sol, graphs[1], varname = :D)
plotgraph(sol, graphs[1], varname = :P)
plotgraph(sol, graphs[1], varname = :Π)













import ModelingToolkit: get_eqs, get_unknowns, get_ps, get_parameter_dependencies, get_observed, get_continuous_events, get_discrete_events, get_defaults, get_systems, get_name, get_iv, get_gui_metadata
get_graphnodes(graph) = PlantModules.nodes(graph)
get_graphnodes(graph::Vector) = reduce(vcat, get_graphnodes.(graph))

function getnodesystem(sol::ODESolution, node)
    system = getfield(sol.prob.f, :sys)
    nodename = string(PlantModules.structmod(node)) * string(PlantModules.id(node))
    nodesystem = [subsys for subsys in getproperty(system, :systems) if get_name(subsys) == Symbol(nodename)][1]
    return nodesystem
end

get_MTKunknown_symbol(s::SymbolicUtils.Symbolic) = s.metadata[ModelingToolkit.VariableSource][2]




indep_values = copy(sol[get_iv(sol.prob.f.sys)]) # values of indepedent variable
append!(indep_values, NaN) # NaN used to cause linebreaks in plot
graphnodes = get_graphnodes(graph)

structmods = PlantModules.structmod.(graphnodes)
if structmod != Symbol("") # user specified e.g. `structmod = :Leaf` => filter out other nodes
    chosen_structmods = structmods .== structmod
    structmods = structmods[chosen_structmods]
    graphnodes = graphnodes[chosen_structmods]
end

nodesystems = getnodesystem.([sol], graphnodes)
varnames = Dict{Symbol, Vector{Symbol}}() # what variables should be plotted per structmod

for structmod_idx in unique(i -> structmods[i], eachindex(structmods))
    _structmod = structmods[structmod_idx]
    varnames[_structmod] = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[structmod_idx])] |> unique

    if varname != Symbol("") # user specified e.g. varname = :W
        varnames[_structmod] = (_varname in varnames[_structmod]) ? [_varname] : Symbol[]
    end
end

varlist = [
    getproperty(nodesystems[nidx], _varname)
    for nidx in eachindex(nodesystems) for _varname in varnames[structmods[nidx]]
]
cumulvarlengths = length.(varlist) |> cumsum

varlocs = Dict{Symbol, Dict{Symbol, Vector{Int64}}}()
nc = 0 # nodecounter
for _structmod in structmods
    if !haskey(varlocs, _structmod)
        varlocs[_structmod] = Dict{Symbol, Vector{Int64}}()
    end

    for _varname in varnames[_structmod]
        structlocs = get!(varlocs[_structmod], _varname, Int64[])

        nc += 1
        varidxs = collect(get(cumulvarlengths, nc-1, 0)+1:cumulvarlengths[nc])
        append!(structlocs, varidxs)
    end
end

varvalues = sol[reduce(vcat, varlist)] |> x -> reduce(hcat, x) |> x -> [x fill(NaN, size(x, 1))]

plots = []

for _varname in unique(vcat(values(varnames)...))
    curr_varlocs = [varlocs[_structmod][_varname] for _structmod in keys(varlocs)] # vector per structmod with indexes of var values
    groups = [fill(_structmod, group_size) for (_structmod, group_size) in zip(keys(varlocs), length.(curr_varlocs))] |> x -> vcat(x...)
    
    ys = varvalues[vcat(curr_varlocs...), :]' |> x -> vcat(x...)
    xs = repeat(indep_values, length(ys) ÷ length(indep_values))

    plot(xs, ys, group = groups, title = "$_varname")
end