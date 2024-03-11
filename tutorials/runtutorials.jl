if split(pwd(), '\\')[end] != "tutorials"
    cd("./tutorials")
end
using Pkg; Pkg.activate(".")
Pkg.develop(path = "..")

using Pluto
Pluto.run()