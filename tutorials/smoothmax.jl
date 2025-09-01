using Plots
using Pkg; Pkg.activate("./tutorials")
using LaTeXStrings

hardmax(x) = max(0, x)
smoothmax(x; k) = x * 1/(1+exp(-k*x))
LSE(x; α) = 1/α * log( sum( [exp(α * i) for i in [x, zero(x)]] ) )

plot_font = "Computer Modern"

begin
    default()
    default(fontfamily=plot_font,
        linewidth=2, framestyle=nothing, label=nothing, grid=true)
    scalefontsizes(1.5)

    plot(title = "Comparison of thresholding functions", xlims = (-0.5, 1), 
        size = (800, 600), xlabel = L"x")
    plot!(hardmax, color = :black, label = L"\mathrm{max}(x, 0)", linewidth = 2)
    plot!(x -> LSE(x, α = 5), color = RGB(1, 0.5, 0.5), linestyle = :dash, label = L"\mathrm{LSE}_5(x, 0)")
    plot!(x -> LSE(x, α = 10), color = RGB(0.75, 0.5, 0.75), linestyle = :dash, label = L"\mathrm{LSE}_{10}(x, 0)")
    plot!(x -> LSE(x, α = 50), color = RGB(0.5, 0.5, 1), linestyle = :dash, label = L"\mathrm{LSE}_{50}(x, 0)")
end

savefig("./tutorials/fig_plantmodules_max.pdf")