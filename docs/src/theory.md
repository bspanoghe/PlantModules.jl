# Theoretical overview

`PlantModules.jl` provides a core functionality for modeling plant growth based on water flows and turgor-driven growth.

## Water flow
Water flow is modelled using the concept of water potentials. We decompose the total water potential $\Psi$ as follows:
```math
\Psi = \Psi_p + \Psi_\pi + \Psi_h + \Psi_m \, ,
```

with the following components:
- The pressure potential $\Psi_p$, equivalent to the turgor pressure $P$.
- The osmotic potential $\Psi_\pi$, equivalent to the negative osmotic pressure $-\Pi$. It expresses the effect of the concentration of osmotically active solutes $M$ in the water and is often approximated using the van 't Hoff equation $\Psi_\pi = - R \, T \, M$, with $R$ the ideal gas constant  \[$8.314$ J mol K\] and $T$ the absolute temperature \[K\].
- The gravitational potential $\Psi_h$, expressing the effects of gravity. For water at a height $h$ \[m\] above a reference level, $\Psi_h = \rho_w \, h \, g$ with $\rho_w$ the density of water \[kg m$^{-3}$\] and $g$ the gravitational acceleration \[N kg$^{-1}$\].
- The matric potential $\Psi_m$, expressing the effects of solid-liquid interactions. It's effect is negligible within plants and is often omitted there, but it plays an important role in the water potential of the soil.

## Turgor-driven growth
Our model for turgor-driven growth is based on the seminal model by Lockart (1965)[^lockhart]
```math
\frac{1}{l}\frac{\mathrm{d}l}{\mathrm{d} t} = \phi \, P \, ,
```
where $l$ is the length of the compartment, $t$ the time, $\phi$ the extensibility and $P$ the turgor pressure.

We make the following extensions:
- Include elastic deformation;
- Include a smooth yield stress threshold;
- Use the equation separately for every axis of the compartment, each with their own mechanical parameters.
This yields the following equation for every axis $i$:
```math
\frac{1}{l_i} \frac{dl_i}{dt} = \phi_i \, (P - P_e)_{\hat{+}_\alpha} + \frac{1}{E_i} \, \frac{\mathrm{d}P}{\mathrm{d}t} \, , 
```
where $P_e$ is the yield stress, $E$ the elastic modulus, and $(x)_{\hat{+}_\alpha} = \mathrm{LSE}_{\alpha}(x, 0)$ with $\mathrm{LSE}_{\alpha}$  the LogSumExp (LSE) function with steepness parameter $\alpha$.

[^lockhart]: Lockhart, J. A. (1965): “An analysis of irreversible plant cell elongation,” Journal of Theoretical Biology, 8 (2), 264–275.