
# ### Photosynthesis

# #### PAR flux from data

PAR_data = readcsv("./tutorials/clouddata/PAR_data.csv")
PAR_interpolation = LinearInterpolation(last.(PAR_data), first.(PAR_data)) # μmol / m^2 / s
plot(PAR_interpolation)

get_PAR_flux(t) = PAR_interpolation(t) / η_photon # J / m^2 / s

# #### carbon assimilation rate

@memoize function get_assimilation_rate(PAR_flux, T, LAI, k)
	Kelvin_to_C = -273.15
	meteo = Atmosphere(T = T + Kelvin_to_C, Wind = 1.0, P = 101.3, Rh = 0.65, Ri_PAR_f = PAR_flux)
	m = ModelList(
		Fvcb(), # calculate CO2 assimilation rate
		Medlyn(0.03, 0.92), # calculate stomatal conductance, see https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2486.2010.02375.x
		Beer(k), # calculate amount of light intercepted
		status = (Tₗ = meteo[:T], LAI = LAI, Cₛ = meteo[:Cₐ], Dₗ = meteo[:VPD], RI_PAR_f = meteo[:Ri_PAR_f])
	)
	run!(m, meteo)
    assimilation_rate = m[:A][1] # extract result of the first (and only) timestep
	return assimilation_rate #! negative values returned
end

@register_symbolic get_assimilation_rate(PAR_flux, T, LAI, k)

plot(t -> get_assimilation_rate(get_PAR_flux(t), 293.15, 8.0, 0.2), xlims = (0.0, 24.0))

# #### module definition
using PlantBiophysics, PlantBiophysics.PlantMeteo, PlantSimEngine
using Memoization

function photosynthesis_module(; name, T, M, shape, M_c)
	@constants (
		uc = (10^-6 * 10^-4 * 60^2), [description = "Unit conversion from (µmol / m^2 / s) to (mol / cm^2 / hr)", unit = u"(mol/cm^2/hr) / (µmol/m^2/s)"],
		    # the output from PlantBiophysics.jl is in different units than we use for our ODEs, so we need to change this
        t_unit = 1, [description = "Dummy constant for correcting units", unit = u"hr"],
	)
	@parameters (
		T = T, [description = "Temperature", unit = u"K"],
		LAI = 8.0, [description = "Leaf Area Index", unit = u"cm^2 / cm^2"],
		k = 0.215, [description = "Light extinction coefficient", unit = u"N/N"],
        M_c = M_c, [description = "Rate of carbon consumption", unit = u"hr^-1"],
	)
	@variables (
        M(t) = M, [description = "Osmotically active metabolite content", unit = u"mol / cm^3"],
		PF(t), [description = "Incoming PAR flux", unit = u"J / s / m^2"],
		A(t), [description = "Carbon assimilation rate", unit = u"µmol / m^2 / s"],
		D(t)[1:getdimensionality(shape)], [description = "Dimensions of compartment", unit = u"cm"],
    )

    eqs = [
		PF ~ get_PAR_flux(t)
		A ~ get_assimilation_rate(PF, T, LAI, k)
        d(M) ~ uc * A * surface_area(shape, D) / volume(shape, D) - M_c*M
    ]
    return System(eqs, t; name, checks = false)
end