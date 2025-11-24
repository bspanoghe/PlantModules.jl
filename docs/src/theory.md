# Theoretical overview

`PlantModules.jl` provides a core functionality for modeling plant growth based on water flows and turgor-driven growth.

## Water flow
Water flow is modelled using the concept of water potentials. We decompose the total water potential $\Psi$ as follows:
$$
\Psi = \Psi_p + \Psi_\pi + \Psi_h + \Psi_m \, ,
$$
with the following components:
- The pressure potential $\Psi_p$, equivalent to the turgor pressure $P$.
- The osmotic potential $\Psi_\pi$, equivalent to the negative osmotic pressure $-\Pi$. It expresses the effect of the concentration of osmotically active solutes $M$ in the water and is often approximated using the van 't Hoff equation $\Psi_\pi = - R \, T \, M$, with $R$ the ideal gas constant  \[$8.314$ J mol K\] and $T$ the absolute temperature \[K\].
- The gravitational potential $\Psi_h$, expressing the effects of gravity. For water at a height $h$ \[m\] above a reference level, $\Psi_h = \rho_w \, h \, g$ with $\rho_w$ the density of water \[kg m$^{-3}$\] and $g$ the gravitational acceleration \[N kg$^{-1}$\].
- The matric potential $\Psi_m$, expressing the effects of solid-liquid interactions. It's effect is negligible within plants and is often omitted there, but it plays an important role in the water potential of the soil.

## Turgor-driven growth
Our model for turgor-driven growth is based on the seminal model by Lockart (1965)[^lockhart]



[^lockhart]: Lockhart, J. A. (1965): “An analysis of irreversible plant cell elongation,” Journal of
Theoretical Biology, 8 (2), 264–275.