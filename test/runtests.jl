module tests
    using Test, PlantModules, PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Plots
    import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv

    @testset "Structure definition" include("./test_plantstructure.jl")
    @testset "System generation" include("./test_generate_system.jl")
    @testset "Plotting" include("./test_plotting.jl")
    @testset "Reading plant structure files" include("test_readgraphs.jl")
end