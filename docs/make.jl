# https://docs.julialang.org/en/v1/manual/documentation/
# https://documenter.juliadocs.org/stable/
# https://adrianhill.de/julia-ml-course/docs/

cd(@__DIR__)

using Pkg; Pkg.activate(".")
using Documenter, PlantModules
# using Literate
using PlutoStaticHTML

bopts = BuildOptions("./src/tutorials", output_format = documenter_output)
build_notebooks(bopts, ["2_real_model_notebook.jl"])

# Literate.markdown("./src/tutorial1.jl", "./src")

pages = [
    "Introduction" => "index.md",
    "Basics tutorial" => "tutorial-1.md",
    "API" => "api.md",
    "Theoretical overview" => "theory.md",
]

makedocs(; sitename = "Whatever", pages, modules = [PlantModules], warnonly = true)

# deploydocs(repo = "https://github.com/bspanoghe/PlantModules.jl")