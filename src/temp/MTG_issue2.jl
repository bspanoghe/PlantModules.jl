using MultiScaleTreeGraph
mtg = Node(MutableNodeMTG("/", "Scene", 0, 0), Dict())
insert_child!(mtg, MutableNodeMTG("/", "Plant", 0, 1))
insert_child!(mtg[1], MutableNodeMTG("/", "Internode", 0, 2))
insert_child!(mtg[1][1], MutableNodeMTG("<", "Leaf", 0, 2))
insert_child!(mtg[1][1], MutableNodeMTG("<", "Leaf", 0, 2))

new_mtg = delete_nodes!(mtg, filter_fun = node -> node_mtg(node).scale != 2)
length(new_mtg) == 3 # false

import MultiScaleTreeGraph: delete_nodes!

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