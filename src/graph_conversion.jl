const MSTG = MultiScaleTreeGraph
const PM = PlantModules

# MultiScaleTreeGraph.jl #

"""
    convert_to_MSTG(graph)

Convert a graph to MultiScaleTreeGraph.jl format. Requires definition of standard graph functions and tree graph functions (see graph_functions.jl)
"""
function convert_to_MSTG(graph)
    rootnode = root(graph)
    MSTGgraph = MSTGify_node(rootnode)
    add_children!(MSTGgraph, rootnode, graph)
    return MSTGgraph
end

function add_children!(MSTGnode::MultiScaleTreeGraph.Node, node::Dict, graph::Dict)
    for chnode in PM.children(node, graph)
        MSTG_chnode = MSTG.Node(MSTGnode, MSTGify_node_MSTG(chnode), MSTGify_node_attributes(chnode))
        add_children!(MSTG_chnode, chnode, graph)
    end
end

MSTGify_node(node::Dict) = MSTG.Node(MutableNodeMSTG("/", string(PM.structmod(node)), 0, 0), PM.attributes(node))
MSTGify_node_MSTG(node::Dict) = MutableNodeMSTG("/", string(PM.structmod(node)), 0, 0)
MSTGify_node_attributes(node::Dict) = PM.attributes(node)

# PlantGraphs.jl #

"""
    convert_to_PG(graph)

Convert a graph to PlantGraphs.jl format. Requires definition of standard graph functions and tree graph functions (see graph_functions.jl)
"""
function convert_to_PG(graph)
    define_nodes!(graph)
    
    rootnode = root(graph)
    PGgraph = PlantGraphs.StaticGraph(PGify_node(rootnode))
	add_children!(PGgraph, PGgraph.insertion, rootnode, graph)
    return PGgraph
end

function add_children!(PGgraph::PlantGraphs.StaticGraph, node_id::Int, node::Dict, graph::Dict)
    for chnode in PlantModules.children(node, graph)
		targeted_append!(PGgraph, PGify_node(chnode), node_id)
		chnode_id = PGgraph.insertion
		add_children!(PGgraph, chnode_id, chnode, graph)
    end
end

PGify_node(node) = eval(PlantModules.structmod(node))(values(PlantModules.attributes(node))...)

function targeted_append!(g::PlantGraphs.StaticGraph, n::PlantGraphs.Node, id::Int)
    nID = PlantGraphs.append!(g, id, PlantGraphs.GraphNode(n))
    PlantGraphs.update_insertion!(g, nID)
    return g
end

function define_nodes!(graph)
	structmods = [
		[PlantModules.structmod(node), collect(keys(PlantModules.attributes(node)))...]
		for node in PlantModules.nodes(graph)
	] |> x -> unique(vec -> vec[1], x)

	for structmod in structmods
		define_node!(structmod...)
	end
end

function define_node!(name, fields...) #! add type annotation to fields?
	expr = quote
		struct placeholdername <: PlantGraphs.Node
			placeholderfield
		end
	end
	expr.args[2].args[2].args[1] = name
	expr.args[2].args[3].args = collect(fields)
	eval(expr)
	return nothing
end