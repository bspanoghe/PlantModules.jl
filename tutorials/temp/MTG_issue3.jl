using MultiScaleTreeGraph
mtg = Node(MutableNodeMTG("/", "foo", 0, 0), Dict())
insert_child!(mtg, MutableNodeMTG("/", "bar", 0, 0))

insert_generation!(mtg, MutableNodeMTG("/", "xyzzy", 0, 0))
parent(mtg[1][1]) == mtg[1]

import MultiScaleTreeGraph: insert_generation!, rechildren!, reparent!, new_node_MTG

function insert_generation!(node::Node{N,A},
    template, attr_fun=node -> A(), maxid=[max_id(node)]) where {N<:AbstractNodeMTG,A}

    maxid[1] += 1

    new_node = Node(
        maxid[1],
        node,
        children(node),
        new_node_MTG(node, template),
        copy(attr_fun(node)),
        Dict{String,Vector{Node{N,A}}}()
    )

    # Add the new node as the only child of the node:
    rechildren!(node, Node{N,A}[new_node])

    # Add the new node as parent of the children
    for chnode in children(new_node)
        setfield!(chnode, :parent, new_node)
    end
    
    return node
end