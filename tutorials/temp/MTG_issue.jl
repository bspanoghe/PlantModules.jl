using MultiScaleTreeGraph
mtg = Node(MutableNodeMTG("/", "foo", 0, 0), Dict())
insert_child!(mtg, MutableNodeMTG("/", "bar", 0, 0))

new_mtg = delete_node!(mtg) # errors

# Move `reparent!` after `link!`
function delete_node_fixed!(node::Node{N,A}; child_link_fun=new_child_link) where {N<:AbstractNodeMTG,A}
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

mtg = Node(MutableNodeMTG("/", "foo", 0, 0), Dict())
insert_child!(mtg, MutableNodeMTG("/", "bar", 0, 0))

new_mtg = delete_node_fixed!(mtg)
new_mtg