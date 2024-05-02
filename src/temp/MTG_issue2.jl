using MultiScaleTreeGraph
mtg = Node(MutableNodeMTG("/", "foo", 0, 0), Dict())
insert_child!(mtg, MutableNodeMTG("/", "bar", 0, 0))
insert_child!(mtg[1], MutableNodeMTG("/", "boo", 0, 0))
insert_child!(mtg[1], MutableNodeMTG("/", "far", 0, 0))
