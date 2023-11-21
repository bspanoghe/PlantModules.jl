#=
Created on 06/10/2023 17:39:25
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

A simple modular plant model using ModelingToolkit

A plant compartment produces sugars according to a first-order kinetics.

It also draws or looses water, depending on its connections
=#

using ModelingToolkit, DifferentialEquations, Plots

@variables t
D = Differential(t)

function plant_compartment(;name)
    @parameters k p r
    @variables t V(t) X(t) C(t) ΔP(t)
    eqs = [
        C ~ X / V,  # concentration is metabolite / volume
        D(X) ~ p * V - r * C * V,  # change in metabolite is photosyn - resp
        D(V) ~ k * ΔP  # change in volume due to turgor presure
    ]
    return ODESystem(eqs; name, defaults = Dict(k => 0.2, p => 1.0, r => 0.3))
end

# make organs
@named root = plant_compartment()
@named branch = plant_compartment()
@named shoot1 = plant_compartment()
@named shoot2 = plant_compartment()

# connect via water potential
connections = [root.ΔP ~ (root.C - 0.2) + (root.C - branch.C),
                branch.ΔP ~ 3branch.C - root.C - shoot1.C - shoot2.C,
                shoot1.ΔP ~ 2shoot1.C - branch.C - 2,
                shoot2.ΔP ~ 2shoot2.C - branch.C - 2]


plant = compose(ODESystem(connections, name = :plant), root, branch, shoot1, shoot2)
plant_simp = structural_simplify(plant)

full_equations(plant_simp)

u0 = [root.X => 0.1,
      root.V => 0.2,
      branch.X => 0.1,
      branch.V => 0.2,
      shoot1.X => 0.1,
      shoot1.V => 0.2,
      shoot2.X => 0.1,
      shoot2.V => 0.2,
      ]
      
prob = ODEProblem(plant_simp, u0, (0.0, 5.0))

plot(solve(prob))