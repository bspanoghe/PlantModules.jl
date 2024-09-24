module basictests
    using Test, PlantModules, PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Plots
    import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv
    @testset "Basic functionality and plotting" include("./generate_system_and_plotting.jl")
end

module readingtests
    using Test, PlantModules
    @testset "Reading plant structure files" include("test_readgraphs.jl")
end