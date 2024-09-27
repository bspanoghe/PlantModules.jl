@independent_variables t, [description = "Time", unit = u"hr"]; #! add documentation
d = Differential(t);

sigm_func(x; A = 1, s = 1, d = 0) = A / (1 + exp(-s*(x - d)))

"""
    hydraulic_module(; name, T, shape, Γ, P, D)

Returns a ModelingToolkit ODESystem describing the turgor-driven growth of a plant compartment.
WARNING: this module still requires an equation to be given for the osmotically active metabolite content M.
"""
function hydraulic_module(; name, T, shape::Shape, Γ, P, D)
    num_D = length(shape.ϵ_D)
    @constants (
        P_0 = 0.0, [description = "Minimum pressure", unit = u"MPa"],
        R = 8.314, [description = "Ideal gas constant", unit = u"MPa * cm^3 / K / mol"], # Pa = J/m^3 => J = Pa * m^3 = MPa * cm^3
        MPa_unit = 1.0, [description = "Dummy constant for correcting units", unit = u"MPa"], #!
    )
    @parameters (
        T = T, [description = "Temperature", unit = u"K"],
        ϵ_D[1:num_D] = shape.ϵ_D, [description = "Dimensional elastic modulus", unit = u"MPa"],
        ϕ_D[1:num_D] = shape.ϕ_D, [description = "Dimensional extensibility", unit = u"MPa^-1 * hr^-1"],
        Γ = Γ, [description = "Critical turgor pressure", unit = u"MPa"],
        ρ_w = 1.0, [description = "Density of water", unit = u"g / cm^3"],
    )
    @variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        Π(t), [description = "Osmotic water potential", unit = u"MPa"],
        P(t) = P, [description = "Hydrostatic potential", unit = u"MPa"],
        M(t), [description = "Osmotically active metabolite content", unit = u"mol / cm^3"], # m^3 so units match in second equation () #! extend validation function so L is ok?
        W(t) = PlantModules.volume(shape, D) * ρ_w, [description = "Water content", unit = u"g"],
        D(t)[1:num_D] = D, [description = "Dimensions of compartment", unit = u"cm"],
        V(t), [description = "Volume of compartment", unit = u"cm^3"],
        ΣF(t), [description = "Net incoming water flux", unit = u"g / hr"],
        
        ΔP(t), [description = "Change in hydrostatic potential", unit = u"MPa / hr"],
        ΔW(t), [description = "Change in water content", unit = u"g / hr"],
        ΔD(t)[1:num_D], [description = "Change in dimensions of compartment", unit = u"cm / hr"],
    )

    eqs = [
        Ψ ~ P - Π, # Water potential consists of a solute- and a pressure component
        Π ~ R*T*M, # Solute component is determined by concentration of dissolved metabolites
        ΔW ~ ΣF, # Water content changes due to flux (depending on water potentials as defined in connections)
        V ~ W / ρ_w, # Volume is directly related to water content  
        V ~ volume(shape, D), # Volume is also directly related to compartment dimensions
        [ΔD[i] ~ D[i] * (ΔP/ϵ_D[i] + ϕ_D[i] * sigm_func((P - Γ)/MPa_unit, A = P - Γ, s = 10)) for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure
        # [ΔD[i] ~ D[i] * (ΔP/ϵ_D[i] + ϕ_D[i] * max(P - Γ, P_0)) for i in eachindex(D)]..., # Compartment dimensions can only change due to a change in pressure

        d(P) ~ ΔP,
        d(W) ~ ΔW,
        [d(D[i]) ~ ΔD[i] for i in eachindex(D)]...,
    ]
    # return ODESystem(eqs, t; name, continuous_events = [P ~ Γ]) #!
    return ODESystem(eqs, t; name)
end

"""
    environmental_module(; name, T, W_max, W_r)

Returns a ModelingToolkit ODESystem describing a non-growing water reservoir.
WARNING: this module still requires an equation to be given for the total water potential Ψ.
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
    return ODESystem(eqs, t; name)
end

"""
    constant_carbon_module(; name, C)

Returns a ModelingToolkit ODESystem describing a constant concentration of osmotically active metabolite content.
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
    return ODESystem(eqs, t; name)
end

"""
    Ψ_soil_module(; name)

Returns a ModelingToolkit ODESystem describing an empirical relationship between
the total water potential of the soil and its relative water content.
"""
function Ψ_soil_module(; name)
	@variables (
        Ψ(t), [description = "Total water potential", unit = u"MPa"],
        W_r(t), [description = "Relative water content", unit = u"g / g"],
    )
    @constants MPa_unit = 1 [description = "Dummy constant for correcting units of empirical equation", unit = u"MPa"]

	eqs = [Ψ ~ MPa_unit * -(1/(100*W_r) + 1) * exp((39.8 - 100*W_r) / 19)]

	return ODESystem(eqs, t; name)
end

"""
    Ψ_air_module(; name, T)

Returns a ModelingToolkit ODESystem describing the relationship between
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

	return ODESystem(eqs, t; name)
end

#! documentation
function sizedep_K_module(; name, K_s, shape::PlantModules.Shape)
    num_D = length(shape.ϵ_D)

    @parameters (
        K_s = K_s, [description = "Specific hydraulic conductivity of compartment", unit = u"g / hr / MPa / cm^2"],
    )
    @variables (
        K(t), [description = "Hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
		D(t)[1:num_D], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		K ~ K_s * PlantModules.cross_area(shape, D)
    ]

    return ODESystem(eqs, t; name)
end

#! documentation
function constant_K_module(; name, K)
    @parameters (
        K_val = K, [description = "Value for the hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
    )
    @variables (
        K(t), [description = "Hydraulic conductivity of compartment", unit = u"g / hr / MPa"],
    )

    eqs = [
		K ~ K_val
    ]

    return ODESystem(eqs, t; name)
end

# Module connections #

## Connection information

"""
    hydraulic_connection(; name)

Returns a ModelingToolkit ODESystem describing a water flow connection between two hydraulics-based functional modules. 
"""
function hydraulic_connection(; name)
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1", unit = u"g / hr"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1", unit = u"g / hr / MPa"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2", unit = u"g / hr / MPa"],
        Ψ_1(t), [description = "Total water potential of compartment 1", unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2", unit = u"MPa"],
    )

    eqs = [
        F ~ min(K_1, K_2) * (Ψ_2 - Ψ_1) #! mean?
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, reverse_order) = [ 
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
		connection_MTK.K_1 ~ node_MTK.K,
        connection_MTK.K_2 ~ nb_node_MTK.K,
    ]
    return ODESystem(eqs, t; name), get_connection_eqset
end

## Node behaviour in function of all their connections

multi_connection_eqs(node_MTK, connection_MTKs) = [
    node_MTK.ΣF ~ sum([connection_MTK.F for connection_MTK in connection_MTKs])
]

# Default values #

default_shape = Sphere(ϵ_D = [1.0], ϕ_D = [0.001])

default_values = Dict(
    hydraulic_module => Dict(:T => 298.15, :shape => default_shape, :Γ => 0.3, :P => 0.5, :D => [15]),
    constant_carbon_module => Dict(:M => 25e-6),
    environmental_module => Dict(:T => 298.15, :W_max => 1e6, :W_r => 1.0),
    Ψ_soil_module => Dict(),
    Ψ_air_module => Dict(:T => 298.15),
    sizedep_K_module => Dict(:K_s => 10, :shape => default_shape),
    constant_K_module => Dict(:K => 5),
    hydraulic_connection => Dict(),
)