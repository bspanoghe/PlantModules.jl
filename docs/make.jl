cd(@__DIR__)

using Pkg; Pkg.activate(".")
using Documenter, DocumenterInterLinks
using PlutoStaticHTML
using PlantModules

rebuild_notebooks = true

if rebuild_notebooks
    cd("./src/tutorials")
    bopts = BuildOptions(".", output_format = documenter_output)
    build_notebooks(bopts)
    cd(@__DIR__)
end

pages = [
    "Introduction" => "index.md",
    "Basics tutorial" => "tutorials/1_basics_notebook.md",
    "Advanced tutorial" => "tutorials/2_advanced_notebook.md",
    "API" => "api.md",
    "List of variables" => "variables.md",
    "Theoretical overview" => "theory.md",
]

links = InterLinks(
    "ModelingToolkit" => "https://docs.sciml.ai/ModelingToolkit/stable/"
)

makedocs(;
    sitename = "PlantModules",
    pages,
    modules = [PlantModules],
    plugins = [links],
    format = Documenter.HTML(size_threshold = 2000 * 1024),
    warnonly = true
)

deploydocs(repo = "github.com/bspanoghe/PlantModules.jl.git")