<<<<<<< HEAD
module tests
    using Test, PlantModules, PlantGraphs, ModelingToolkit, OrdinaryDiffEq, Plots
    import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv

    @testset "Structure definition" include("./test_plantstructure.jl")
    @testset "System generation" include("./test_generate_system.jl")
    @testset "Plotting" include("./test_plotting.jl")
    @testset "Reading plant structure files" include("test_readgraphs.jl")
=======
module basictests
    using Test, PlantModules, PlantGraphs, ModelingToolkit, DifferentialEquations, Plots
    import ModelingToolkit: get_eqs, get_systems, get_unknowns, get_defaults, get_name, get_iv
    @testset "Basic functionality and plotting" include("./generate_system_and_plotting.jl")
end

module readingtests
    using Test, PlantModules
    @testset "Reading plant structure files" include("./test_readgraphs.jl")
>>>>>>> c71c803cbf5e84f69a5edcf0ed7a38b1cb5d4b7b
end