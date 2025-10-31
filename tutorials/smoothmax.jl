using Plots, LaTeXStrings

hardmax(x) = max(0, x)
smoothmax(x; k) = x * 1/(1+exp(-k*x))
LSE(x; α) = 1/α * log( sum( [exp(α * i) for i in [x, zero(x)]] ) )

# plot_font = "Computer Modern"
include(homedir() * raw"\Documents\Github\Caverns_of_code\Julia\Lifehacks\catpuccin\get_palette.jl")
cpalette = get_palette("prettycolors")
plotdir = homedir() * raw"\Documents\Github\Den_of_evil\Non-note files\images\\"

begin
    default()
    default(linewidth=2)
    scalefontsizes(1.5)

    plot(title = "Comparison of thresholding functions", xlims = (-0.5, 1), 
        size = (800, 600), xlabel = L"x")
    plot!(hardmax, color = :black, label = L"\mathrm{max}(x, 0)", linewidth = 2)
    plot!(x -> LSE(x, α = 5), color = cpalette[1], linestyle = :dash, label = L"\mathrm{LSE}_5(x, 0)")
    plot!(x -> LSE(x, α = 10), color = cpalette[2], linestyle = :dash, label = L"\mathrm{LSE}_{10}(x, 0)")
    plot!(x -> LSE(x, α = 50), color = cpalette[3], linestyle = :dash, label = L"\mathrm{LSE}_{50}(x, 0)")
end

savefig(plotdir * "fig_plantmodules_max.pdf")