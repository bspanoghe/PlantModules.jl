# For functions I don't need anymore (but what if I end up needing them later ?!) #

# Apply a function to a node and all its descendants
function iteratedescendants(node, graph, func::Function; kwargs...)
	func(node; kwargs...)
	for chnode in PlantModules.children(node, graph)
		iteratedescendants(chnode, graph, func; kwargs...)
	end
end

# Default behaviour: start from graph root
iteratedescendants(graph, func::Function; kwargs...) = iteratedescendants(PlantModules.root(graph), graph, func; kwargs...)

# Fill in MTK_systems, MTK_connections, MTK_u0 for the given node
function add_MTK_info!(node; model_defaults, module_defaults, module_coupling, struct_connections,
	func_connections, MTK_systems, MTK_connections, MTK_u0)

	MTKsystem = getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params)
	# MTKconnections = getMTKconnections(node) # needs all connected nodes to be defined as MTKsystems already
end

"""
	root(graph)

Returns the root node of the given graph.
"""
root(graph) = error("Function not yet defined for given input type.")


"""
	children(node)

Returns the children nodes of the given node.
"""
children(node) = error("Function not yet defined for given input type.")

root(graph::PlantGraphs.StaticGraph) = graph[graph.root]
root(graph::PlantGraphs.GraphNode) = graph

children(node::PlantGraphs.GraphNode, graph::PlantGraphs.StaticGraph) = [graph[child_id] for child_id in node.children_id]
0
## iteratedescendants 
testgraph = Foo(1) + (Foo(2), Foo(3) + (Foo(7), Foo(9)))
bars = Int[]
iteratedescendants(testgraph, (x; extra) -> push!(bars, x.data.bar + extra), extra = 3)
bars == [1, 2, 3, 7, 9] .+ 3


# get MTK systems and the vector of equations defining connections between MTK systems of all nodes between two graphs as specified in the intergraph connections
function get_intergraph_connection_info(graphs, intergraph_connection, func_connections,
	get_connection_eqs, default_params, default_u0s, MTK_system_dicts)

	graph1, graph2 = graphs[intergraph_connection[1]]
	MTK_system_dict, nb_MTK_system_dict = MTK_system_dicts[intergraph_connection[1]]

	nodes = [node for node in PlantModules.nodes(graph1) if PlantModules.nodetype(node) in intergraph_connection[2]]
	nb_nodes = [node for node in PlantModules.nodes(graph2) if PlantModules.nodetype(node) in intergraph_connection[2]]
	
	if isempty(nodes) || isempty(nb_nodes) # no neighbours no connection info
		return ODESystem[], Equation[]
	end

	full_connection_MTKs = Vector{ODESystem}(undef, length(nodes) * length(nb_nodes)) # we're fully connecting the valid nodes from both graphs
	full_connection_equations = Equation[]

	for (node_nr, node) in enumerate(nodes)
		node_MTK = MTK_system_dict[PlantModules.id(node)]
		
		nb_node_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
		connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes))

		for (nb_nr, nb_node) in enumerate(nb_nodes)
			nb_node_MTKs[nb_nr] = nb_MTK_system_dict[PlantModules.id(nb_node)]
			connection_MTKs[nb_nr] = get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
		end
	
		push!(full_connection_equations, get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs)...)
		full_connection_MTKs[node_nr:node_nr+length(nb_nodes)-1] = connection_MTKs
	end

	return full_connection_MTKs, full_connection_equations
end

intergraph_connections = struct_connections[2]
for intergraph_connection in intergraph_connections
	intergraph_connection_MTKs, intergraph_connection_equations = get_intergraph_connection_info(graphs, intergraph_connection, func_connections,
		get_connection_eqs, default_params, default_u0s, MTK_system_dicts
	)
	push!(connection_MTKs, intergraph_connection_MTKs...)
	push!(connection_equations, intergraph_connection_equations...)
end


function get_symmetry_info(connection_dict, get_symmetry_eqs)
	symmetry_eqs = Equation[]
	for connection_key in keys(connection_dict)
		rev_key = reverse(connection_key) # key of reverse connection (by construction of keys)
		append!(symmetry_eqs, get_symmetry_eqs(connection_dict[connection_key], connection_dict[rev_key]))
	end

	return symmetry_eqs
end

get_symmetry_eqs(connection_MTK, reverse_connection_MTK) = [
	connection_MTK.F ~ -reverse_connection_MTK.F
] #! put this in func_modules and explain to user they gotta provide this stuff mr white yo

append!(connection_equations, get_symmetry_info(connection_dict, get_symmetry_eqs))

# get MTK systems and the vector of equations defining connections between MTK systems of a node and its neighbour nodes
function get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
	node_MTK = MTK_system_dicts[graphnr][PlantModules.id(node)]

	nb_node_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
	connection_MTKs = Vector{ODESystem}(undef, length(nb_nodes))
	connection_keys = Vector{Tuple{Tuple, Tuple}}(undef, length(nb_nodes))

	for (nb_nr, (nb_node, nb_node_graphnr)) in enumerate(zip(nb_nodes, nb_node_graphnrs))
		nb_node_MTKs[nb_nr] = MTK_system_dicts[nb_node_graphnr][PlantModules.id(nb_node)]
		connection_MTKs[nb_nr] = get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
		connection_keys[nb_nr] = ((graphnr, PlantModules.id(node)), (nb_node_graphnr, PlantModules.id(nb_node)))
	end

	connection_equations = get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs)
	connection_dict = Dict(Pair.(connection_keys, connection_MTKs))
	return connection_dict, connection_equations
end





# Functional modules #

@variables t, [description = "Time", unit = u"hr"]; #! add documentation
d = Differential(t);

"""
    hydraulic_module(; name, T, ρ_w, shape, Γ)

Returns a ModelingToolkit ODESystem describing the turgor-driven growth of a plant compartment.
WARNING: this module still requires an equation to be given for the osmotically active metabolite content M.
"""
function hydraulic_module(; name, T, ρ_w, shape::Shape, Γ, P, W, D)
    num_D = length(shape.ϵ_D)
    @constants (
        P_0 = 0.0, [description = "Minimum pressure", unit = u"MPa"],
        R = 8.314e-6, [description = "Ideal gas constant", unit = u"MJ / K / mol"]
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / m^3"],
        ϵ_D[1:num_D] = shape.ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        ϕ_D[1:num_D] = shape.ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        Γ = Γ, [description = "Critical turgor pressure", unit = u"MPa"]
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / m^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
        W(t) = W, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"m"],
        V(t), [description = "Shape of compartment", unit = u"m^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"m / hr"],
    )

    eqs = [
        Ψ ~ Π + P, # Water potential consists of a solute- and a pressure component
        Π ~ -R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        V ~ W / ρ_w, # Shape is directly related to water content        
        V ~ volume(shape, D), # Shape is also directly related to compartment dimensions
        [ΔD[i] ~ D[i] * (ΔP/ϵ_D[i] + ϕ_D[i] * LSE(P - Γ, P_0, γ = 100)) for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return ODESystem(eqs, t; name)
end

"""
    constant_carbon_module(; name, C)

Returns a ModelingToolkit ODESystem describing constant osmotically active metabolite content.
"""
function constant_carbon_module(; name, M_const)
    @parameters M_const = M_const [description = "Value of constant M concentration", unit = u"mol / m^3"]
    @variables M(t) [description = "Osmotically active metabolite content", unit = u"mol / m^3"]

    eqs = [
        M ~ M_const,
    ]
    return ODESystem(eqs, t; name)
end

"""
    environmental_module(; name, T, ρ_w, W_max)

Returns a ModelingToolkit ODESystem describing a non-growing water reservoir.
WARNING: this module still requires an equation to be given for the total water potential Ψ.
"""
function environmental_module(; name, T, ρ_w, W_max, W_r)
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ρ_w = ρ_w, [description = "Density of water", unit = u"g / m^3"],
        W_max = W_max, [description = "Water capacity of compartment", unit = u"g"],
        )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        W(t) = W_r * W_max, [description = "Water content", unit = u"g"],
        W_r(t), [description = "Relative water content", unit = u"g / g"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],

        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
    )

    eqs = [
        W_r ~ W / W_max,
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        d(W) ~ ΔW,

        # prevent variables from being erased from existence
        Ψ ~ Ψ,
        T ~ T
    ]
    return ODESystem(eqs, t; name)
end

# Module connections #

"""
    hydraulic_connection(; name, K)

Returns a ModelingToolkit ODESystem describing a water flow connection between two hydraulics-based functional modules. 
"""
function hydraulic_connection(; name, K)
    @parameters (
        K = K, [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1)
    ]
    return ODESystem(eqs, t; name)
end

# Helper functions #

## Unitful is a dangerous beast
val(x) = x
val(x::Quantity) = x.val

## Forbidden rites #! try fixing this... creative solution?
import Base.exp
import Base.log
exp(x::Quantity) = exp(val(x))*unit(x)
log(x::Quantity) = log(val(x))*unit(x)

"""
    LSE(x, y, ...; γ = 1)

LogSumExp, a smooth approximation for the maximum function. 
The temperature parameter γ determines how close the approximation is to the actual maximum.
This implementation makes use of the log-sum-exp trick (https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations) to ensure numerical stability.
"""
LSE(x::Real...; γ = 1) = log(sum(exp.(γ .* x .- maximum(x))) ) / γ + maximum(x)
@register_symbolic LSE(x)

# Default values #

default_params = (
    hydraulic_module = (
        T = 298.15, ρ_w = 1.0e6, shape = Sphere(ϵ_D = [1.0], ϕ_D = [1.0]), Γ = 0.3
    ),
    constant_carbon_module = (
        M_const = 0.2,
    ),
    environmental_module = (
        T = 298.15, ρ_w = 1.0e6, W_max = 1e9
    ),
    hydraulic_connection = (
        K = 600,
    )
)

default_u0s = (
    hydraulic_module = (
        P = 0.1, M = 200.0, D = [0.1], W = volume(Sphere(ϵ_D = [1.0], ϕ_D = [1.0]), [0.1]) * 1.0e6,
    ),
    constant_carbon_module = (
    ),
    environmental_module = (
        W_r = 1,
    ), 
    hydraulic_connection = (
    )
)










function constant_carbon_module(; name, M, shape::Shape, D)
    @parameters ( #! change input to M and calculate corresponding M_amount?
        M_amount = M * PlantModules.volume(shape, D), [description = "Amount of osmotically active metabolite content", unit = u"mol"],
    )

    @variables (
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / m^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
        V(t), [description = "Volume of compartment", unit = u"m^3"],
    )

    eqs = [
        M ~ M_amount/V
    ]
    return ODESystem(eqs, t; name)
end


mutable struct Root <: Node 
	D::Vector
end

mutable struct PNode <: Node
end
mutable struct Internode <: Node
	D::Vector
end
mutable struct Meristem <: Node
	G
end
mutable struct Leaf <: Node
	D::Vector
end

struct Soil <: Node end
struct Air <: Node end

## parameters
branching_prob = 0.1
G_th = 5
G0 = 10
rootsize0 = [0.05, 0.5]
internodesize0 = [0.05, 0.5]
leafsize0 = [0.5, 0.2, 0.01]

## rules

shoot_rule = Rule(Meristem, lhs = meristem -> rand() < data(meristem).G, rhs = meristem -> PNode() + (Meristem(branching_prob) + Leaf(leafsize0), Internode(internodesize0) + Meristem(1)))
# root_rule = Rule(Meristem, lhs = meristem -> rand() < data(meristem).G, rhs = meristem -> Root(rootsize0) + Meristem(branching_prob), Root(rootsize0) + Meristem(branching_prob), Root(rootsize0) + Meristem(branching_prob))

## grow 'em
axiom = PNode() + Internode(internodesize0) + Meristem(G0)

plant_graph = Graph(axiom = axiom, rules = (shoot_rule,))
[rewrite!(plant_graph) for _ in 1:10]












function plotgraph(sol::ODESolution, graph; struct_module::Symbol = Symbol(""), func_varname::Symbol = Symbol(""))
    if struct_module == Symbol("") && func_varname == Symbol("") # Ya get nothin'
        struct_modules = PlantModules.nodes(graph) .|> PlantModules.structmod |> unique
        nodesystems_by_structmod = [getnodesystems(sol, graph, struct_module) for struct_module in struct_modules]
        
        func_varnames = [
            [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[1])]
            for nodesystems in nodesystems_by_structmod
        ]  |> x -> reduce(vcat, x) |> unique
        nodesystems = reduce(vcat, nodesystems_by_structmod)

        plots = []

        for func_varname in func_varnames
            myplot = Plots.plot();
            for nodesystem in nodesystems #! oh the inefficiency
                if func_varname in [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystem)]
                    getplot!(sol, nodesystem, func_varname);
                end
            end
            push!(plots, myplot)
        end
        return plots

    elseif func_varname == Symbol("") # Only structural module provided
        nodesystems = getnodesystems(sol, graph, struct_module)

        func_varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[1])] |> unique
    
        plots = []
        for func_varname in func_varnames
            myplot = Plots.plot();
            for nodesystem in nodesystems
                getplot!(sol, nodesystem, func_varname);
            end
            push!(plots, myplot)
        end
        return plots
    elseif struct_module == Symbol("") # Only variable name provided
        struct_modules = PlantModules.nodes(graph) .|> PlantModules.structmod |> unique
        nodesystems = [getnodesystems(sol, graph, struct_module) for struct_module in struct_modules] |> x -> reduce(vcat, x)
        myplot = Plots.plot();
        for nodesystem in nodesystems
            if func_varname in get_MTKunknown_symbol.(get_unknowns(nodesystem))
                getplot!(sol, nodesystem, func_varname);
            end
        end
        return(myplot)
    else # Everything provided
        nodesystems = getnodesystems(sol, graph, struct_module)
        myplot = Plots.plot();
        for nodesystem in nodesystems
            if func_varname in get_MTKunknown_symbol.(get_unknowns(nodesystem))
                getplot!(sol, nodesystem, func_varname);
            end
        end
        return(myplot)
    end
end

function plotgraph(sol::ODESolution, graphs::Vector; struct_module::Symbol = Symbol(""), func_varname::Symbol = Symbol(""))
    if struct_module == Symbol("") && func_varname == Symbol("") # Ya get nothin'
        nodesystems_by_structmod = []
        for graph in graphs
            struct_modules = PlantModules.nodes(graph) .|> PlantModules.structmod |> unique
            append!(nodesystems_by_structmod, [getnodesystems(sol, graph, struct_module) for struct_module in struct_modules])
        end

        func_varnames = [[get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[1])]
            for nodesystems in nodesystems_by_structmod
        ]  |> x -> reduce(vcat, x) |> unique
        nodesystems = reduce(vcat, nodesystems_by_structmod)
        plots = []

        for func_varname in func_varnames
            myplot = Plots.plot();
            for nodesystem in nodesystems #! oh the inefficiency
                if func_varname in [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystem)]
                    getplot!(sol, nodesystem, func_varname);
                end
            end
            push!(plots, myplot)
        end
        return plots

    elseif func_varname == Symbol("") # Only structural module provided
        nodesystems = [getnodesystems(sol, graph, struct_module) for graph in graphs] |> x -> reduce(vcat, x)

        func_varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystems[1])] |> unique
    
        plots = []
        for func_varname in func_varnames
            myplot = Plots.plot();
            for nodesystem in nodesystems
                getplot!(sol, nodesystem, func_varname);
            end
            push!(plots, myplot)
        end
        return plots
    elseif struct_module == Symbol("") # Only variable name provided
        nodesystems = []
        for graph in graphs
            struct_modules = PlantModules.nodes(graph) .|> PlantModules.structmod |> unique
            append!(nodesystems, [getnodesystems(sol, graph, struct_module) for struct_module in struct_modules] |> x -> reduce(vcat, x))
        end

        myplot = Plots.plot();
        for nodesystem in nodesystems
            if func_varname in get_MTKunknown_symbol.(get_unknowns(nodesystem))
                getplot!(sol, nodesystem, func_varname);
            end
        end
        return(myplot)
    else # Everything provided
        nodesystems = [getnodesystems(sol, graph, struct_module) for graph in graphs] |> x -> reduce(vcat, x)
        myplot = Plots.plot();
        for nodesystem in nodesystems
            if func_varname in get_MTKunknown_symbol.(get_unknowns(nodesystem))
                getplot!(sol, nodesystem, func_varname);
            end
        end
        return(myplot)
    end
end









import MultiScaleTreeGraph: delete_nodes!, delete_nodes_!

function delete_node!(node::Node{N,A}; child_link_fun=new_child_link) where {N<:AbstractNodeMTG,A}
    if isroot(node)
        if length(children(node)) == 1
            # If it has only one child, make it the new root:
            chnode = children(node)[1]
            # Add to the new root the mandatory root attributes:
            root_attrs = Dict(
                :symbols => node[:symbols],
                :scales => node[:scales],
                :description => node[:description]
            )

            append!(chnode, root_attrs)

            link!(chnode, child_link_fun(chnode))
            reparent!(chnode, nothing)

            node_return = chnode
        else
            error("Can't delete the root node if it has several children")
        end
    else
        parent_node = parent(node)

        if !isleaf(node)
            # We re-parent the children to the parent of the node.
            for chnode in children(node)
                # Updating the link of the children:
                link!(chnode, child_link_fun(chnode))
                addchild!(parent_node, chnode; force=true)
            end
        end

        # Delete the node as child of his parent:
        deleteat!(children(parent_node), findfirst(x -> node_id(x) == node_id(node), children(parent_node)))
        node_return = parent_node
    end

    node = nothing

    return node_return
end

function delete_nodes!(
    node;
    scale=nothing,
    symbol=nothing,
    link=nothing,
    all::Bool=true, # like continue in the R package, but actually the opposite
    filter_fun=nothing,
    child_link_fun=new_child_link
	)

    # Check the filters once, and then compute the descendants recursively using `descendants_`
    check_filters(node, scale=scale, symbol=symbol, link=link)
    filtered = is_filtered(node, scale, symbol, link, filter_fun)

    while filtered
        node = delete_node!(node)
        filtered = is_filtered(node, scale, symbol, link, filter_fun)

        # Don't go further if all == false
        !all && return
    end

    delete_nodes!_(node, scale, symbol, link, all, filter_fun, child_link_fun)

    return node
end

function delete_nodes!_(node, scale, symbol, link, all, filter_fun, child_link_fun)
    if !isleaf(node)
        # First we apply the algorithm recursively on the children:
        chnodes = children(node)
        nchildren = length(chnodes)
        #? Note: we don't use `for chnode in chnodes` because it may delete dynamically during traversal, so we forget to traverse some nodes
        for chnode in chnodes[1:nchildren]
            delete_nodes!_(chnode, scale, symbol, link, all, filter_fun, child_link_fun)
        end
    end

    # Then we work on the node itself. This ensures that its children will not be deleted
    # afterwards (the deletion is acropetal, i.e. from leaves to root)

    # Is there any filter happening for the current node? (true is deleted):
    filtered = is_filtered(node, scale, symbol, link, filter_fun)

    if filtered
        delete_node!(node, child_link_fun=child_link_fun)
    end
end



















# based on https://link.springer.com/article/10.1007/BF00195075
irradiance_data = [max(0, 2000 * sin(t/24*2*pi - 8) + randn()) for t in 0:10*24]
@memoize get_irradiance(t) = irradiance_data[floor(Int, t+1)] + (t-floor(t)) * irradiance_data[ceil(Int, t+1)]
# get_irradiance(t) = max(0, 2000 * sin(t/24*2*pi - 8))
@register_symbolic get_irradiance(t)

#=
P(I; ϕ, P_m, Θ) = (ϕ*I + P_m - ((ϕ*I + P_m)^2 - 4*Θ*ϕ*I*P_m)^0.5) / (2*Θ)
eV = 1.602176634e-19 # joules per electron
N_A = 6.02214076e23 # particles per mole
μmol_e_to_J = 10^-6 * eV * N_A
P_J(PF; ϕ, P_m, Θ) = P(μmol_e_to_J*PF; ϕ, P_m, Θ)
xl = 1000
plot(x -> P_J(x; ϕ = 0.9, P_m = 33, Θ = 0.8), xlims = (0, xl));
scatter!(0:xl, get_assimilation_rate.(0:xl, 293.15, 1, 0.7))
=#


function photosynthesis_module(; name, M, shape)
	@constants (
		uc1 = (10^-6 * 60^2), [description = "Unit conversion from (µmol/s) to (mol/hr)", unit = u"(mol/hr) / (µmol/s)"],
	)
	@parameters (
		ϕ = 0.9, [description = "Maximum quantum yield", unit = u"mol/mol"],
		P_m = 33, [description = "Light-saturated photosynthetic rate", unit = u"µmol / m^2 / s"],
		Θ = 0.8, [description = "Convexity of photosynthetic rate"],
		carbon_decay_rate = 0.1, [description = "Rate at which carbon is consumed for growth", unit = u"hr^-1"],
	)

	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / m^3"], # m^3 so units match in second equation (Pa = J/m^3) #! extend validation function so L is ok?
		I(t), [description = "Irradiance", unit = u"μmol / s / m^2"], #! make sure not to use variable name from other func mod used in same struct mod
		PR(t), [description = "Photosynthetic rate", unit = u"µmol / m^2 / s"], # = carbon assimilation rate
		D(t)[1:length(shape.ϵ_D)], [description = "Dimensions of compartment", unit = u"m"],
    )

    eqs = [
		I ~ get_irradiance(t)
		PR ~ (ϕ*I + P_m - ((ϕ*I + P_m)^2 - 4*Θ*ϕ*I*P_m)^0.5) / (2*Θ)
        d(M) ~ uc1 * leafarea(shape, D) * PR / volume(shape, D) - carbon_decay_rate*M # convert µmol => mol and s^-1 => hr^-1
    ]
    return ODESystem(eqs, t; name, checks = false) #! checks back to true?
end

# n_samples = 1000
# lower_bound = [0.0, 273.15]
# upper_bound = [600.0, 303.15]
# PAR_Ts = sample(n_samples, lower_bound, upper_bound, RandomSample())
# As = get_assimilation_rate.(first.(PAR_Ts), last.(PAR_Ts), 2.0, 0.5)
# surr = SecondOrderPolynomialSurrogate(PAR_Ts, As, lower_bound, upper_bound)

# test_PAR_Ts = sample(n_samples, lower_bound, upper_bound, RandomSample())
# test_As = get_assimilation_rate.(first.(test_PAR_Ts), last.(test_PAR_Ts), 2.0, 0.5)
# mean((test_As - surr.(test_PAR_Ts)).^2)
# @btime surr([400.2, 295.3])

# get_est_assimilation_rate(PAR, T) = surr([PAR, T])
# @register_symbolic get_est_assimilation_rate(PAR, T)






# Helper functions #

## Unitful is a dangerous beast
val(x) = x
val(x::Quantity) = x.val

## Forbidden rites #! try fixing this... creative solution?
import Base.exp
import Base.log
exp(x::Quantity) = exp(val(x))*unit(x)
log(x::Quantity) = log(val(x))*unit(x)

"""
    LSE(x, y, ...; γ = 1)

LogSumExp, a smooth approximation for the maximum function. 
The temperature parameter γ determines how close the approximation is to the actual maximum.
"""
LSE(x::Real...; γ = 1) = log(sum(exp.(γ .* x )) ) / γ
@register_symbolic LSE(x)