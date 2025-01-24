using Pkg; Pkg.activate("./tutorials")
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqDefault

function my_module(; name, W_r, W_max)
    @parameters W_max = W_max
    @variables W(t) = W_r * W_max

    eqs = [
        D(W) ~ -W / W_max,
    ]
    return ODESystem(eqs, t; name)
end

sys = my_module(; :name => :sys, :W_max => 10, :W_r => 0.5)
defaults(sys)
sys = my_module(; :name => :sys, :W_max => 5, :W_r => 0.5)
defaults(sys)