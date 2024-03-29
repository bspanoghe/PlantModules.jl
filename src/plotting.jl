#! Add Plots kwargs and maybe add RecipesBase.jl (get rid of Plots dependency)?

# Get the symbol representation of a MTK unknown (variable)
get_MTKunknown_symbol(s::SymbolicUtils.Symbolic) = s.metadata[ModelingToolkit.VariableSource][2]

function getnodesystem(sol::ODESolution, node)
    system = getfield(sol.prob.f, :sys)
    nodename = string(PlantModules.structmod(node)) * string(PlantModules.id(node))
    nodesystem = [subsys for subsys in getproperty(system, :systems) if get_name(subsys) == Symbol(nodename)][1]
    return nodesystem
end

function getnodesystems(sol::ODESolution, graph, struct_module::Symbol)
    matching_nodes = [node for node in PlantModules.nodes(graph) if PlantModules.structmod(node) == struct_module]
    if isempty(matching_nodes)
        error("No structural module $(struct_module) found in graph")
    end
    return [getnodesystem(sol, node) for node in matching_nodes]
end

function getplot(sol::ODESolution, nodesystem, func_varname::Symbol)
    func_var = getproperty(nodesystem, func_varname)
    myplot = getplot(sol, func_var)
    return myplot
end

function getplot(sol::ODESolution, func_var::Num)
    myplot = Plots.plot(sol, idxs = func_var);
    return myplot
end

function getplot(sol::ODESolution, func_var::Symbolics.Arr)
    myplot = Plots.plot();
    for sub_func_var in func_var
        Plots.plot!(sol, idxs = sub_func_var, label = "$sub_func_var");
    end
    return myplot
end

function getplot!(sol::ODESolution, nodesystem, func_varname::Symbol)
    func_var = getproperty(nodesystem, func_varname)
    getplot!(sol, func_var)
end

function getplot!(sol::ODESolution, func_var::Num)
    Plots.plot!(sol, idxs = func_var);
end

function getplot!(sol::ODESolution, func_var::Symbolics.Arr)
    for sub_func_var in func_var
        Plots.plot!(sol, idxs = sub_func_var, label = "$sub_func_var");
    end
end

"""
    plotnode(sol::ODESolution, node; func_varname::Symbol)

Returns a plot for every functional variable of the given node for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable.
"""
function plotnode(sol::ODESolution, node; func_varname::Symbol = Symbol(""))
    nodesystem = getnodesystem(sol, node)

    if func_varname == Symbol("") # No functional variable name given => show all
        func_varnames = [get_MTKunknown_symbol(unknown) for unknown in get_unknowns(nodesystem)] |> unique

        plots = [getplot(sol, nodesystem, func_varname) for func_varname in func_varnames]
        return plots
    else
        return getplot(sol, nodesystem, func_varname)
    end
end

"""
    plotgraph(sol::ODESolution, graph; struct_module::Symbol, func_varname::Symbol)

Returns a plot for every functional variable for every node of the given graph for the given solution `sol`.
Optionally, the user can give the name of a functional variable to only return a plot of this variable,
give the name of a structural module to limit considered nodes to those of this type, or both.
"""
function plotgraph(sol::ODESolution, graph; struct_module::Symbol = Symbol(""), func_varname::Symbol = Symbol(""))
    if struct_module == Symbol("") && func_varname == Symbol("") # Ya get nothin'
        struct_modules = PlantModules.nodes(graph) .|> PlantModules.structmod |> unique
        nodesystems_by_structmod = [getnodesystems(sol, graph, struct_module) for struct_module in struct_modules]
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