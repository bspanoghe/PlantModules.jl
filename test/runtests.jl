module tests1
    using Test, PlantModules, PlantGraphs, ModelingToolkit
    import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv
    @testset "Basic functionality" include("./basictests.jl") 
end

module tests2
    using Test, PlantModules, PlantGraphs, ModelingToolkit, DifferentialEquations, Plots
    import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv
    @testset "Basic functionality and plotting" include("./lotka_volterra_tests.jl")
end