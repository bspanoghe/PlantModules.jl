"""
    readXEG(file::String)

Reads in the graph from a given XEG format file.
"""
function readXEG(file::String)
    lines = readlines(file)
    graph = Dict{Int, Dict}()

    curr_id = 0
    root_parent_id = -1
    parents = Dict{Int, Int}()

    for line in lines
        line_element = split(line)[1]
        if line_element == "<node"
            id = htmlmatch("id", line) |> x -> parse(Int, x)
            structmod = htmlmatch("type", line)

            graph[id] = Dict(:id => id, :parent_id => 0, :child_ids => Int[], :attributes => Dict{Symbol, Any}(), :type => Symbol(structmod))
            curr_id = id
        elseif line_element == "<edge"
            id = htmlmatch("src_id", line) |> x -> parse(Int, x)
            child_id = htmlmatch("dest_id", line) |> x -> parse(Int, x)
            push!(graph[id][:child_ids], child_id)

            parents[child_id] = id
        elseif line_element == "<property" && occursin("value=", line)
            name = htmlmatch("name", line)
            val = htmlmatch("value", line)
            if !isnothing(match(r"^[0-9]+$", val))
                val = parse(Int, val)
            elseif !isnothing(match(r"^-?[0-9]*\.[0-9]+(E-[0-9]+)?$", val))
                val = parse(Float64, val)
            end
            graph[curr_id][:attributes][Symbol(name)] = val
        end
    end

    for node in graph
        node[2][:parent_id] = get(parents, node[2][:id], root_parent_id)
    end

    return graph
end

htmlmatch(element, line) = match(Regex("(?<=$(element)=\").+?(?=\")"), line).match