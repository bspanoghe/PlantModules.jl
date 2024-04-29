struct MyPGNode{T} <: PlantGraphs.Node
    attributes::Dict
end

MyPGNode(type::Symbol, attributes) = MyPGNode{type}(attributes)