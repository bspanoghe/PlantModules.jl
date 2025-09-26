# # Setup
@independent_variables t, [description = "Time", unit = u"hr"]; # independent variable
d = Differential(t); # differential operator

# # Modules

# ## Hydraulics
"""
    hydraulic_module(; name, shape, ϕ_D, ϵ_D, Γ, T, D, Ψ, M)

Return a ModelingToolkit System describing the turgor-driven growth of a plant compartment.

This module still requires a module describing the osmotically active metabolite content M.
"""
function hydraulic_module(; name, shape::Shape, ϕ_D, ϵ_D, Γ, T, D, Ψ, M)
    D, ϕ_D, ϵ_D = [correctdimensionality(shape, var) for var in [D, ϕ_D, ϵ_D]] 
        # turns scalar values into vectors of correct length

    num_D = getdimensionality(shape)
    R = 8.314
    P = Ψ + M*R*T

    @constants (
        R = R, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
        P_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"],
        ρ_w = 1.0, [description = "Density of water", unit = u"g / cm^3"],
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϕ_D[1:num_D] = ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        ϵ_D[1:num_D] = ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        Γ = Γ, [description = "Yield turgor pressure", unit = u"MPa"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation ()
        W(t) = volume(shape, D) / ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr", guess = 0.0], #!
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P + Π, # Water potential consists of a solute- and a pressure component
        Π ~ -R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        V ~ W / ρ_w, # Volume is directly related to water content  
        V ~ volume(shape, D), # Volume is also directly related to compartment dimensions
        [ΔD[i] ~ D[i]*ϕ_D[i]*P_unit*logsumexp((P - Γ)/P_unit, α = 100) + D[i]*ΔP/ϵ_D[i] for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    return System(eqs, t; name)
end

"""
    environmental_module(; name, T, W_max, W_r)

Return a ModelingToolkit System describing a non-growing water reservoir.

This module still requires a module describing the total water potential Ψ.
"""
function environmental_module(; name, T, W_max, W_r)
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        W_max = W_max, [description = "Water capacity of compartment", unit = u"g"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        W(t) = W_r * W_max, [description = "Water content", unit = u"g"],
        W_r(t), [description = "Relative water content", unit = u"g / g"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],

        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
    )

    eqs = [
        W_r ~ W / W_max,
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        d(W) ~ ΔW,
    ]
    return System(eqs, t; name)
end

# ## Carbon dynamics

"""
    simple_photosynthesis_module(; name, M, shape)

Return a ModelingToolkit System describing a concentration of osmotically active metabolite content increasing during day and decreasing in function of the current concentration.
"""
function simple_photosynthesis_module(; name, M, shape, t_sunrise, t_sunset, A_max, M_c)
    @constants (
        t_unit = 1, [description = "Dummy constant for correcting units", unit = u"hr"],
    )
    @parameters (
        t_sunrise = t_sunrise, [description = "Time of sunrise (hours past midnight)", unit = u"hr"],
        t_sunset = t_sunset, [description = "Time of sunset (hours past midnight)", unit = u"hr"],
        A_max = A_max, [description = "Maximum carbon assimilation rate", unit = u"mol / cm^2 / hr"],
        M_c = M_c, [description = "Rate of carbon consumption", unit = u"hr^-1"],
	)
	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / cm^3"],
		A(t), [description = "Carbon assimilation rate", unit = u"mol / cm^2 / hr"],
		D(t)[1:getdimensionality(shape)], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		A ~ smooth_daynight(t/t_unit, t_sunrise/t_unit, t_sunset/t_unit, zero(A_max), A_max, smoothing = 1.0)
        d(M) ~ A * surface_area(shape, D) / volume(shape, D) - M_c * M
    ]
    return System(eqs, t; name)
end

"""
    constant_carbon_module(; name, C)

Return a ModelingToolkit System describing a constant concentration of osmotically active metabolite content.
"""
function constant_carbon_module(; name, M)
    @parameters (
        M_value = M, [description = "Constant value of osmotically active metabolite content", unit = u"mol / cm^3"],
    )

    @variables (
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"],
    )

    eqs = [
        M ~ M_value
    ]
    return System(eqs, t; name)
end

# ## Water potentials (for parts that don't use `hydraulic_module`)

"""
    Ψ_air_module(; name, T)

Return a ModelingToolkit System describing the relationship between
the total water potential of the air and its relative water content.
"""
function Ψ_air_module(; name, T)
	@variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        W_r(t), [description = "Relative water content", unit = u"g / g"],
    )
	@parameters T = T [description = "Temperature", unit = u"K"]
	@constants (
        R = 8.314, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"],
        V_w = 18, [description = "Molar volume of water", unit = u"cm^3/mol"]
    )

	eqs = [Ψ ~ R * T / V_w * log(W_r)] # Spanner equation (see e.g. https://academic.oup.com/insilicoplants/article/4/1/diab038/6510844)

	return System(eqs, t; name)
end

"""
    Ψ_soil_module(; name)

Return a ModelingToolkit System describing an empirical relationship between
the total water potential of the soil and its relative water content.
"""
function Ψ_soil_module(; name)
	@variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        W_r(t), [description = "Relative water content", unit = u"g / g"],
    )
    @constants Ψ_unit = 1 [description = "Dummy constant for correcting units", unit = u"MPa"]

	eqs = [Ψ ~ Ψ_unit * soilfunc(W_r)]

	return System(eqs, t; name)
end

soilfunc(W_r) = -(1/W_r) * exp(-30*W_r) # empirical equation for soil water potential


# ## Hydraulic conductivity

"""
    K_module(; name, K_s, shape::Shape)

Return a ModelingToolkit System describing the hydraulic conductance of a 
compartment as the product of its specific hydraulic conductance and an area of the compartment.
"""
function K_module(; name, K_s, shape::Shape, area_func::Function)
    num_D = getdimensionality(shape)

    @parameters (
        K_s = K_s, [description = "Specific hydraulic conductivity", unit = u"g / hr / MPa / cm^2"],
    )
    @variables (
        K(t), [description = "Hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
		D(t)[1:num_D], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		K ~ K_s * area_func(shape, D)
    ]

    return System(eqs, t; name)
end

"""
    constant_K_module(; name, K_s, shape::Shape)

Return a ModelingToolkit System describing the hydraulic conductance of a 
compartment as a constant.
"""
function constant_K_module(; name, K)
    @parameters (
        K_value = K, [description = "Hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
    )
    @variables (
        K(t), [description = "Hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
    )

    eqs = [
		K ~ K_value
    ]

    return System(eqs, t; name)
end

# # Module connections

"""
    hydraulic_connection(; name)

Returns a ModelingToolkit System describing a water flow connection between two hydraulics-based compartments.

This module assumes the compartments have a specified hydraulic conductivities.
"""
function hydraulic_connection(; name)
    @constants (
        K_unit = 1, [description = "Dummy constant for correcting units", unit = u"g / hr / MPa"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        K(t), [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1", unit = u"g / hr / MPa"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2", unit = u"g / hr / MPa"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1),
        K ~ min(K_1, K_2)
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = [ 
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
		connection_MTK.K_1 ~ node_MTK.K,
        connection_MTK.K_2 ~ nb_node_MTK.K,
    ]
    return System(eqs, t; name), get_connection_eqset
end

"""
    evaporation_connection(; name, t_sunrise, t_sunset, η_night)

Returns a ModelingToolkit System describing a water flow connection between two hydraulics-based compartments that decreases at night.

This module assumes the compartments have a specified hydraulic conductivities.
"""
function evaporation_connection(; name, t_sunrise, t_sunset, η_night)
    @constants (
        t_unit = 1, [description = "Dummy constant for correcting units", unit = u"hr"],
        K_unit = 1, [description = "Dummy constant for correcting units", unit = u"g / hr / MPa"],
    )
    @parameters (
        t_sunrise = t_sunrise, [description = "Time of t_sunrise (hours past midnight)", unit = u"hr"],
        t_sunset = t_sunset, [description = "Time of t_sunset (hours past midnight)", unit = u"hr"],
        η_night = η_night, [description = "Relative hydraulic conductivity at night", unit = u"(g / hr / MPa) / (g / hr / MPa)"]
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        K(t), [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1", unit = u"g / hr / MPa"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2", unit = u"g / hr / MPa"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1) * smooth_daynight(t/t_unit, t_sunrise/t_unit, t_sunset/t_unit, η_night),
        K ~ min(K_1, K_2)
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = 
        [ 
            connection_MTK.Ψ_1 ~ node_MTK.Ψ,
            connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
            connection_MTK.K_1 ~ node_MTK.K,
            connection_MTK.K_2 ~ nb_node_MTK.K,
        ]
    
    return System(eqs, t; name), get_connection_eqset
end

"""
    constant_hydraulic_connection(; name, K)

Returns a ModelingToolkit System describing a water flow connection between two hydraulics-based functional modules.

This module specifies a constant hydraulic conductivity between the compartments.
"""
function constant_hydraulic_connection(; name, K)
    @parameters (
        K = K, [description = "Hydraulic conductivity of connection", unit = u"g / hr / MPa"],
    )

    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ K * (Ψ_2 - Ψ_1)
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = [ 
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
    ]
    return System(eqs, t; name), get_connection_eqset
end

# # Node behaviour in function of all their connections

multi_connection_eqs(node_MTK, connection_MTKs) = [
    node_MTK.ΣF ~ sum([connection_MTK.F for connection_MTK in connection_MTKs])
        # "node's net incoming water flux equals the sum of all connected water flows"
]

# # Default values

default_values = Dict(
    :shape => Cylinder(), :ϕ_D => 0.02, :ϵ_D => 50.0, :Γ => 0.3,
    :T => 298.15, :D => [0.5, 5.0], :Ψ => 0.0, :M => 300e-6, :W_max => 1e6, :W_r => 0.5, 
    :K_s => 10.0, :K => 1e3, :t_sunrise => 8, :t_sunset => 20, :η_night => 0.1, :A_max => 2e-6, 
    :M_c => 0.05, :area_func => cross_area
)