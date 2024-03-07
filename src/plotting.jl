function getnodesystem(sol::ODESolution, node)
    system = sol.prob.f.sys
    nodename = string(PlantModules.nodetype(node)) * string(PlantModules.id(node))
    nodesystem = [subsys for subsys in system.systems if subsys.name == Symbol(nodename)][1]
    return nodesystem
end

function getnodesystems(sol::ODESolution, graph, struct_module::Symbol)
    matching_nodes = [node for node in PlantModules.nodes(graph) if PlantModules.nodetype(node) == struct_module]
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

function plot(sol::ODESolution; node) #! Give other name?
    nodesystem = getnodesystem(sol, node)

    func_varnames = [unknown.val.metadata[ModelingToolkit.VariableSource][2] for unknown in nodesystem.unknowns] |> unique

    plots = [getplot(sol, nodesystem, func_varname) for func_varname in func_varnames]
    for myplot in plots
        display(Plots.plot(myplot))
        sleep(0.1)
    end
end

function plot(sol::ODESolution; node, func_varname::Symbol)
    nodesystem = getnodesystem(sol, node)
    display(getplot(sol, nodesystem, func_varname))
end

function plot(sol::ODESolution; graph, struct_module::Symbol, func_varname::Symbol)
    nodesystems = getnodesystems(sol, graph, struct_module)
    myplot = Plots.plot();
    for nodesystem in nodesystems
        getplot!(sol, nodesystem, func_varname);
    end
    display(myplot)
end