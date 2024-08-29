using Pkg; Pkg.activate("./tutorials")
include("../../src/PlantModules.jl"); using .PlantModules
using PlantGraphs
using ModelingToolkit, DifferentialEquations, Unitful
using Plots

import ModelingToolkit: t
d = Differential(t)

struct Transect <: Node
    T
end
struct Border <: Node
    T
end

make_transect(T_ins) = Transect(T_ins[1]) + Transect(T_ins[2]) + Transect(T_ins[3]) +
    Transect(T_ins[4]) + Transect(T_ins[5]) + Transect(T_ins[6]) + Transect(T_ins[7]) +
    Transect(T_ins[8]) + Transect(T_ins[9]) + Transect(T_ins[10])

transect = make_transect(20:0.1:20.9)
border = Border(28)
graphs = [transect, border]
intergraph_connections = [[1, 2] => (PlantModules.nodes(transect)[1], :Border)] #! select nodes
struct_connections = [graphs, intergraph_connections]


outside_temp(t) = 4*cos((t - 15)/12*pi) + 23
plot(outside_temp, xlims = (0, 48), xticks = 0:6:48)

daytime_activator(t) = (t%24 >= 8) && (t%24 <= 20)
# plot(daytime_activator, xlims = (0, 48))
@register_symbolic daytime_activator(t)

input_radiation(t) = 400*cos((t - 14)/12*pi) * daytime_activator(t)
plot(input_radiation, xlims = (0, 48))
@register_symbolic input_radiation(t)

function forest_module(; name, T, outside_conductance, cooling_factor)
    @variables T(t) = T T_out(t) ΣF(t)
    @parameters(
        outside_conductance = outside_conductance,
        cooling_factor = cooling_factor,
        )

    eqs = [
        d(T) ~ outside_conductance * (T_out - T) + ΣF - cooling_factor * daytime_activator(t) + 0.01*input_radiation(t),
        T_out ~ outside_temp(t)
    ]

    return ODESystem(eqs, t; name, checks = false)
end

function border_module(; name)
    @variables T(t) ΣF(t)

    eqs = [T ~ outside_temp(t), ΣF ~ ΣF]

    return ODESystem(eqs, t; name, checks = false)
end

function heat_connection(; name, transect_conductance)
    @parameters transect_conductance = transect_conductance
    @variables F T1 T2

    eqs = [
        F ~ transect_conductance * (T2 - T1)
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order) = [ 
        connection_MTK.T1 ~ node_MTK.T,
        connection_MTK.T2 ~ nb_node_MTK.T,
    ]

    return ODESystem(eqs, t; name, checks = false), get_connection_eqset
end

module_coupling = [
    forest_module => [:Transect]
    border_module => [:Border]
]

connecting_modules = [
    (:Border, :Transect) => (heat_connection, []),
    (:Transect, :Transect) => (heat_connection, []),
]

func_connections = [connecting_modules, PlantModules.multi_connection_eqs]

default_params = (
    forest_module = (outside_conductance = 0.5, cooling_factor = 2),
    border_module = (),
    heat_connection = (transect_conductance = 1,)
)

default_u0s = (
    forest_module = (T = 0,),
    border_module = (),
    heat_connection = ()
)

module_defaults = (
    Transect = (),
    Border = (),
)

system = PlantModules.generate_system(default_params, default_u0s,
	module_defaults, module_coupling, struct_connections, func_connections, checkunits = false
)

sys_simpl = structural_simplify(system);
prob = ODEProblem(sys_simpl, ModelingToolkit.missing_variable_defaults(sys_simpl), (0.0, 2*24))
@time sol = solve(prob);

begin
    plot(sol, color = :red);
    plot!(outside_temp, label = "Outside temperature")
end