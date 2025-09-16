const MTG = MultiScaleTreeGraph
const PM = PlantModules

# MultiScaleTreeGraph.jl #

"""
    convert_to_MTG(graph)

Convert a graph to MultiScaleTreeGraph.jl format. Requires definition of standard graph functions and tree graph functions (see graph_functions.jl)
"""
function convert_to_MTG(graph)
    rootnode = getroot(graph)
    MTGgraph = MTGify_node(rootnode)
    add_children!(MTGgraph, rootnode, graph)
    return MTGgraph
end

function add_children!(MTGnode::MultiScaleTreeGraph.Node, node, graph)
    for chnode in PM.getchildren(node, graph)
        MTG_chnode = MTG.Node(MTGnode, MTGify_node_MTG(chnode), MTGify_node_attributes(chnode))
        add_children!(MTG_chnode, chnode, graph)
    end
end

MTGify_node(node) = MTG.Node(MutableNodeMTG("<", string(PM.getstructmod(node)), 0, 0), PM.getattributes(node))
MTGify_node_MTG(node) = MutableNodeMTG("<", string(PM.getstructmod(node)), 0, 0)
MTGify_node_attributes(node) = PM.getattributes(node)

# PlantGraphs.jl #

"""
    convert_to_PG(graph)

Convert a graph to PlantGraphs.jl format. Requires definition of standard graph functions and tree graph functions (see graph_functions.jl)
"""
function convert_to_PG(graph)
    rootnode = PlantModules.getroot(graph)
    PGgraph = PlantGraphs.StaticGraph(PGify_node(rootnode))
	add_children!(PGgraph, PGgraph.insertion, rootnode, graph)
    return PGgraph
end

function add_children!(PGgraph::PlantGraphs.StaticGraph, node_id::Int, node, graph)
    for chnode in PlantModules.getchildren(node, graph)
		targeted_append!(PGgraph, PGify_node(chnode), node_id)
		chnode_id = PGgraph.insertion
		add_children!(PGgraph, chnode_id, chnode, graph)
    end
end

PGify_node(node) = MyPGNode(PlantModules.getstructmod(node), PlantModules.getattributes(node))

function targeted_append!(g::PlantGraphs.StaticGraph, n::PlantGraphs.Node, id::Int)
    nID = PlantGraphs.append!(g, id, PlantGraphs.GraphNode(n))
    PlantGraphs.update_insertion!(g, nID)
    return g
end