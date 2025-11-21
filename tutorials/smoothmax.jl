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

    plot(title = "Comparison of thresholding functions", xlims = (-0.5, 0.5), 
        size = (800, 600), xlabel = L"x")
    plot!(hardmax, color = :black, label = L"\mathrm{max}(x, 0)", linewidth = 2)

    αs = [4, 40, 400]
    for (i, α) in enumerate(αs)
        plot!(x -> LSE(x; α), color = cpalette[i], linestyle = :dash, label = L"\mathrm{LSE}_{%$(α)}(x, 0)")
    end
end
plot!()

# savefig(plotdir * "fig_plantmodules_max.pdf")