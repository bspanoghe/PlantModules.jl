cd(@__DIR__)

using Pkg; Pkg.activate(".")
using Documenter, Literate, PlantModules
Literate.markdown("./src/tutorial1.jl", "./src")

pages = [
    "Introduction" => "index.md",
    "Basics tutorial" => "tutorial1.md",
    "API" => "api.md",
    "Theoretical overview" => "theory.md",
]

makedocs(; sitename = "Whatever", pages, modules = [PlantModules], warnonly = true)

# deploydocs(repo = "https://github.com/bspanoghe/PlantModules.jl")