const MTG = MultiScaleTreeGraph
const PM = PlantModules

MTGify_node(node::Dict) = MTG.Node(MutableNodeMTG("/", string(PM.structmod(node)), 0, 0), PM.attributes(node))
MTGify_node_MTG(node::Dict) = MutableNodeMTG("/", string(PM.structmod(node)), 0, 0)
MTGify_node_attributes(node::Dict) = PM.attributes(node)

function add_children!(MTGnode::MultiScaleTreeGraph.Node, node::Dict, graph::Dict)
    for childnode in PM.children(node, graph)
        MTG_childnode = MTG.Node(MTGnode, MTGify_node_MTG(childnode), MTGify_node_attributes(childnode))
        add_children!(MTG_childnode, childnode, graph)
    end
end

function convert_to_MTG(graph::Dict)
    rootnode = [node for node in PM.nodes(graph) if node[:parent_id] == -1][1]
    MTGgraph = MTGify_node(rootnode)
    add_children!(MTGgraph, rootnode, graph)
    return MTGgraph
end