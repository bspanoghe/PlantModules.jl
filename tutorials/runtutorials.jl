if split(pwd(), '\\')[end] != "tutorials"
    cd("./tutorials")
end
using Pkg; Pkg.activate(".")
using Pluto
Pluto.run()