#=
Created on 21/11/2023 15:31:34
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Checking out with ModelingToolkit how to correctly use dimensions to model volumetric change
=#

using ModelingToolkit, Plots, DifferentialEquations

# soft version of the maximum
σ(x; κ=2) = log(exp(κ * x) + 1)

@variables t

D = Differential(t)

# A TOMATO


@mtkmodel Sphere begin
    @parameters begin 
        ϵᵣ = 8.0
        ϕᵣ = 0.3 
        Γ = 1.0
    end
    @variables begin
        V(t)
        P(t)
        r(t)
    end
    @equations begin
        V ~ 4π/3 *  r^3  # hoe zorgen hier geen beginwaarden?
        D(r) / r ~ D(P) / ϵᵣ + σ(P-Γ)  # deze geef je mee voor elke richting
        D(P) ~  1 - 0.1P  # te bepalen
    end
end

@named tomato = Sphere(;r=0.1, P=0.3, V=4π/3*0.1^3)

equations(tomato)

tomato_simp = structural_simplify(tomato)

sol = ODEProblem(tomato_simp, [D(tomato_simp.r)=>0], (0.0, 5.0))  |> solve
plot(sol)

# CELL

@mtkmodel Cube begin
    @parameters begin 
        ϵ = 8.0
        ϕ = 0.3 
        Γ = 1.0
    end
    @variables begin
        V(t)
        P(t)
        x(t)
        y(t)
        z(t)
    end
    @equations begin
        V ~ x * y * z  # hoe zorgen hier geen beginwaarden?
        D(x) / x ~ D(P) / ϵ + σ(P-Γ)  # deze geef je mee voor elke richting
        D(y) / y ~ D(P) / ϵ + σ(P-Γ)  # deze geef je mee voor elke richting
        D(z) / z ~ D(P) / ϵ + σ(P-Γ)  # deze geef je mee voor elke richting
        D(P) ~  1 - 0.1P  # te bepalen
    end
end