# https://docs.julialang.org/en/v1/manual/documentation/
# https://documenter.juliadocs.org/stable/
# https://adrianhill.de/julia-ml-course/docs/

cd(@__DIR__)

using Pkg; Pkg.activate(".")
using Documenter, PlantModules
using PlutoStaticHTML

cd("./src/tutorials")
bopts = BuildOptions(".", output_format = documenter_output)
build_notebooks(bopts)
cd(@__DIR__)

pages = [
    "Introduction" => "index.md",
    "Basics tutorial" => "tutorials/1_basics_notebook.md",
    "Advanced tutorial" => "tutorials/2_advanced_notebook.md",
    "API" => "api.md",
    "Theoretical overview" => "theory.md",
]

makedocs(; 
    sitename = "PlantModules", 
    pages, 
    modules = [PlantModules],
    format = Documenter.HTML(size_threshold = 2000 * 1024),
    warnonly = true
)

# deploydocs(repo = "https://github.com/bspanoghe/PlantModules.jl")