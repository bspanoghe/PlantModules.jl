using Pkg; Pkg.activate(".")
# include("../src/PlantModules.jl")
using Test, PlantModules, PlantGraphs, ModelingToolkit, DifferentialEquations, Plots

@testset "Set 1" begin
    mutable struct Foo <: Node
        bar::Int
    end

    mutable struct Oof <: Node
        rab::String
    end

    function qux(; name, bar, baz, quux, corge, xyzzy)
        @parameters bar = bar baz = baz quux = quux corge = corge
        @variables t xyzzy(t) = xyzzy
        eqs = [
            bar + baz ~ quux * xyzzy + corge
        ]
        return ODESystem(eqs, t; name)
    end

    function snoo(; name, snee)
        @parameters snee = snee
        @variables t xyzzy1(t) xyzzy2(t) snaw(t)
        eqs = [
            snaw ~ xyzzy1 * xyzzy2 + snee
        ]
        return ODESystem(eqs, t; name)
    end

    PlantModules.attributes(node::Foo) = Dict([:bar => node.bar])
    PlantModules.structmod(::Foo) = :Foo
    PlantModules.structmod(::Oof) = :Oof
    PlantModules.id(::Foo) = 1

    func_module = qux
    structmodule = :Foo
    func_connections = [() => snoo]
    module_coupling = [qux => [:Foo]]
    default_params = (qux = (bar = 100, baz = 100, quux = 100, corge = 100), snoo = (snee = 100,))
    default_u0s = (qux = (xyzzy = 100,),)
    model_defaults = (baz = 42, quux = 42)
    module_defaults = (Foo = (baz = 10, xyzzy = 10),)


    ## getnodeparams

    node = Foo(1)

    @test PlantModules.getnodeparams(node, structmodule, func_module, module_defaults, model_defaults, default_params) ==
        Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100])

    ## getnodeu0s

    @test PlantModules.getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s) ==
        Dict([:xyzzy => 10])

    ## getMTKsystem

    node = Foo(1)
    testsystem = PlantModules.getMTKsystem(node, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
    @test testsystem.name == :(Foo1)
    @test sort(Symbol.(keys(testsystem.defaults))) == sort(collect(keys(Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100, Symbol("xyzzy(t)") => 10]))))
    @test sort(collect(values(testsystem.defaults))) == sort(collect(values(Dict([:bar => 1, :baz => 10, :quux => 42, :corge => 100, Symbol("xyzzy(t)") => 10]))))

    ## get_func_connection
    node = Foo(1)
    nb_node = Foo(2)
    testsystem = PlantModules.get_func_connection(node, nb_node, func_connections, default_params, default_u0s)
    @test testsystem.name == :Foo1_Foo1
    @test sort(Symbol.(keys(testsystem.defaults))) == sort(collect(keys(Dict([:snee => 100]))))
    @test sort(collect(values(testsystem.defaults))) == sort(collect(values(Dict([:snee => 100]))))

    ## get_connection_info
    graphnr = 1
    graph = Foo(123) + (Foo(22), Foo(300))
    node = PlantModules.nodes(graph)[1]
    nb_nodes = PlantModules.nodes(graph)[2:3]
    nb_node_graphnrs = [1, 1]
    MTK_system_dicts = [Dict([PlantModules.nodes(graph)[nodenr].self_id =>
        qux(name = Symbol("Foo$nodenr"), bar = nodenr, baz = nodenr, quux = nodenr, corge = nodenr, xyzzy = nodenr)
        for nodenr in 1:3]
    )]
    get_connection_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
        [connection_MTK.xyzzy1 ~ node_MTK.xyzzy for connection_MTK in connection_MTKs]..., 
        [connection_MTK.xyzzy2 ~ nb_node_MTKs.xyzzy for (connection_MTK, nb_node_MTKs) in zip(connection_MTKs, nb_node_MTKs)]..., 
    ]

    connection_MTKs, connection_equations = PlantModules.get_connection_info(node, graphnr, nb_nodes,
        nb_node_graphnrs, func_connections, get_connection_eqs, default_params, default_u0s, MTK_system_dicts
    )

    @test connection_MTKs[1] == snoo(; name = Symbol("Foo" * string(PlantModules.id(node)) * "_" *
    "Foo" * string(PlantModules.id(nb_nodes[1]))), snee = 100)
    @test connection_MTKs[2] == snoo(; name = Symbol("Foo" * string(PlantModules.id(node)) * "_" *
    "Foo" * string(PlantModules.id(nb_nodes[2]))), snee = 100)

    @test connection_equations == get_connection_eqs(MTK_system_dicts[1][PlantModules.id(node)],
        [MTK_system_dicts[1][PlantModules.id(nb_nodes[1])], MTK_system_dicts[1][PlantModules.id(nb_nodes[2])]],
        connection_MTKs
    )

    ## get_intergraph_neighbours
    graphs = [Oof("1") + (Foo(2), Foo(3), Foo(4)), Oof("!")]
    node_graphnr = 2
    node = PlantModules.nodes(graphs[node_graphnr])[1]

    ### 1
    intergraph_connections = [[1, 2] => (:Foo, :Oof)]
    ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
    @test ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
    @test ig_neighbour_graphnrs == [1, 1, 1]

    ### 2
    intergraph_connections = [[2, 1] => (:Foo, :Oof)]
    ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
    @test isempty(ig_neighbours)
    @test isempty(ig_neighbour_graphnrs)

    ### 3
    intergraph_connections = [[2, 1] => (:Oof, :Foo)]
    ig_neighbours, ig_neighbour_graphnrs = PlantModules.get_intergraph_neighbours(node, node_graphnr, graphs, intergraph_connections)
    @test ig_neighbours == PlantModules.nodes(graphs[1])[2:4]
    @test ig_neighbour_graphnrs == [1, 1, 1]

    ## get_nb_nodes
    graphs = [Oof("1") + (Foo(2), Foo(3), Foo(4)), Foo(90)]
    graphnr = 1
    node = PlantModules.nodes(graphs[graphnr])[1]
    intergraph_connections = [[1, 2] => (:Oof, :Foo)]

    nb_nodes = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
    @test nb_nodes[1] == vcat(PlantModules.nodes(graphs[1])[2:4], PlantModules.nodes(graphs[2])[1])
    @test nb_nodes[2] == [1, 1, 1, 2]
end



@testset "Lotka-Volterra tests" begin
    
    function lotka_volterra(; name, α, β, δ, γ, N, P)
        @parameters α = α β = β δ = δ γ = γ
        @variables t N(t) = N P(t) = P ΣF_N(t) ΣF_P(t)
        d = Differential(t);
        eqs = [
            d(N) ~ α*N - β*N*P + ΣF_N,
            d(P) ~ δ*N*P - γ*P + ΣF_P
        ]
        return ODESystem(eqs, t; name)
    end

    function fountain_of_rabbits(; name, η, N, P)
        @variables t N(t) = N P(t) = P ΣF_N(t) ΣF_P(t)
        d = Differential(t);
        eqs = [
            d(N) ~ η + ΣF_N,
            d(P) ~ ΣF_P
        ]
        return ODESystem(eqs, t; name)
    end

    function wandering_animals(; name, κ)
        @parameters (
            κ = κ,
        )
        @variables (
            t,
            F_N(t),
            N_1(t),
            N_2(t),
            F_P(t),
            P_1(t),
            P_2(t)
        )

        eqs = [
            F_N ~ κ * (N_2 - N_1),
            F_P ~ κ * (P_2 - P_1)
        ]
        return ODESystem(eqs, t; name)
    end

    wandering_eqs(node_MTK, nb_node_MTKs, connection_MTKs) = [
        [connection_MTK.N_1 ~ node_MTK.N for connection_MTK in connection_MTKs]...,
        [connection_MTK.P_1 ~ node_MTK.P for connection_MTK in connection_MTKs]...,
        [connection_MTK.N_2 ~ nb_node_MTK.N for (connection_MTK, nb_node_MTK) in zip(connection_MTKs, nb_node_MTKs)]...,
        [connection_MTK.P_2 ~ nb_node_MTK.P for (connection_MTK, nb_node_MTK) in zip(connection_MTKs, nb_node_MTKs)]...,
        node_MTK.ΣF_N ~ sum([connection_MTK.F_N for connection_MTK in connection_MTKs]),
        node_MTK.ΣF_P ~ sum([connection_MTK.F_P for connection_MTK in connection_MTKs])
    ]

    Base.@kwdef mutable struct Forest <: Node
        N::Int = 20
        P::Int = 10
        δ::Float64 = 1.0
    end

    Base.@kwdef mutable struct Grassland <: Node
        N::Int = 100
        P::Int = 5
    end

    Base.@kwdef mutable struct Cave <: Node end

    default_params = (
        lotka_volterra = (α = 1.5, β = 2.0, δ = 1.5, γ = 0.8),
        fountain_of_rabbits = (η = 10,),
        wandering_animals = (κ = 0.01,)
    )
    default_u0s = (
        lotka_volterra = (N = 30, P = 10),
        fountain_of_rabbits = (N = 1, P = 0),
        wandering_animals = ()
    )
    model_defaults = (β = 1.9)
    module_defaults = (
        Grassland = (β = 0.1,),
        Forest = (δ = 2.3,),
        Cave = (β = 0,)
    )
    module_coupling = [
        lotka_volterra => [:Grassland, :Forest],
        fountain_of_rabbits => [:Cave]
    ]
    graph = Grassland() + Grassland(N = 50, P = 3) + (Forest(δ = 1.8), Forest(N = 5, P = 0, δ = 2.3))
    graph2 = Cave() + Forest() + Grassland()
    graphs = [graph, graph2]
    struct_connections = [graphs, [[1, 2] => (:Forest, :Cave)]] #! make vector of graphs as input unnecessary when theres only 1 graph
    func_connections = [[() => wandering_animals], wandering_eqs]

    graphs, intergraph_connections = struct_connections
    connecting_modules, get_connection_eqs = func_connections


    # graphfuncs #

    node1, node2, node3, node4 = collect(values(graph.nodes))

    ## nodes
    @test issetequal(PlantModules.nodes(graph), [node1, node2, node3, node4])

    ## neighbours
    @test issetequal(PlantModules.neighbours(node2, graph), [node1, node3, node4])

    ## attributes
    @test issetequal(PlantModules.attributes(node4), [:P => 0, :N => 5, :δ => 2.3])

    ## structmod
    @test PlantModules.structmod(node1) == :Grassland

    ## id
    @test allunique(PlantModules.id.([node1, node2, node3, node4]))


    # generate_system #

    ## getMTKsystem
    node1, node2, node3, node4 = collect(values(graph.nodes))

    sys1 = PlantModules.getMTKsystem(node1, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
    @test sys1.name == Symbol(string(PlantModules.structmod(node1)) * string(PlantModules.id(node1)))

    ## get_MTK_system_dicts
    MTK_system_dicts = PlantModules.get_MTK_system_dicts(graphs, module_coupling, module_defaults, model_defaults, default_params, default_u0s)
    @test length(MTK_system_dicts) == 2
    @test length(MTK_system_dicts[1]) == 4
    @test length(MTK_system_dicts[2]) == 3

    ## get_nb_nodes
    node = node4
    graphnr = 1
    nb_nodes, nb_node_graphnrs = PlantModules.get_nb_nodes(node, graphs, graphnr, intergraph_connections)
    @test issetequal(nb_nodes, [node2, PlantModules.nodes(graphs[2])[1]])

    ## get_connection_info
    node = node4
    graphnr = 1
    PlantModules.get_connection_info(node, graphnr, nb_nodes, nb_node_graphnrs, connecting_modules, get_connection_eqs, default_params, default_u0s, MTK_system_dicts)
    #! test

    ## getnodeu0s
    node = node4
    structmodule = :Forest
    func_module = lotka_volterra
    nodeu0s = PlantModules.getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)
    @test issetequal(nodeu0s, [:P => 0, :N => 5])

    node = PlantModules.nodes(graph2)[1]
    structmodule = :Cave
    func_module = fountain_of_rabbits
    nodeu0s = PlantModules.getnodeu0s(node, structmodule, func_module, module_defaults, model_defaults, default_u0s)
    @test issetequal(nodeu0s, [:P => 0, :N => 1])

    ## generate_system
    sys = PlantModules.generate_system(model_defaults, module_defaults, module_coupling, struct_connections, func_connections,
        default_params = default_params, default_u0s = default_u0s
    )

    sys_simpl = structural_simplify(sys);

    prob = ODEProblem(sys_simpl, [], (0, 48))
    sol = solve(prob)

    # plot functions #
    @test PlantModules.plotnode(sol, PlantModules.nodes(graphs[2])[1]) isa Vector
    @test PlantModules.plotnode(sol, PlantModules.nodes(graphs[1])[1], func_varname = :N) isa Plots.Plot{T} where {T}

    @test PlantModules.plotgraph(sol, graphs[1]) isa Vector
    @test PlantModules.plotgraph(sol, graphs[1], struct_module = :Forest) isa Vector
    @test PlantModules.plotgraph(sol, graphs[1], func_varname = :ΣF_P) isa Plots.Plot{T} where {T}
    @test PlantModules.plotgraph(sol, graphs[1], struct_module = :Grassland, func_varname = :P) isa Plots.Plot{T} where {T}

    @test PlantModules.plotgraph(sol, graphs) isa Vector
    @test PlantModules.plotgraph(sol, graphs, struct_module = :Forest) isa Vector
    @test PlantModules.plotgraph(sol, graphs, func_varname = :ΣF_P) isa Plots.Plot{T} where {T}
    @test PlantModules.plotgraph(sol, graphs, struct_module = :Grassland, func_varname = :P) isa Plots.Plot{T} where {T}
end