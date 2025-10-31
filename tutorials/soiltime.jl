using Plots
soil_data = readlines("./tutorials/clouddata/soil_data.csv") |> x -> split.(x, ", ") .|> x -> parse.(Float64, x)
Ψs, Wrs = (first.(soil_data), last.(soil_data) ./ 100)

soilfunc(W_r; a, b) = -(a/W_r) * exp(-b*W_r)

a_est = 3.5
b_est = 5.5

ss_res = sum((Ψs - soilfunc.(Wrs, a = a_est, b = b_est)).^2)
ss_tot = sum((Ψs .- Ψs/length(Ψs)).^2)
r2 = 1 - ss_res/ss_tot

scatter(Wrs, Ψs)
plot!(x -> soilfunc(x, a = 3.5, b = 5.5))