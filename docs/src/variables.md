The following variables are used in the pre-implemented functionality used by the package.

| Name                  | Units                                                | Description                                 | Default value | Reference       |
| --------------------- | ---------------------------------------------------- | ------------------------------------------- | ------------- | --------------- |
| t                     | h                                                    | Time                                        | N/A           | N/A             |
| $\Psi$                | MPa                                                  | Total water potential                       | 0.0           | DeSchepper[^52] |
| $\Pi$                 | MPa                                                  | Osmotic water potential                     | N/A           | N/A             |
| $P$                   | MPa                                                  | Hydrostatic potential                       | N/A           | N/A             |
| M                     | mol cm$^{-3}$                                        | Osmotically active metabolite content       | 300e-6        | None            |
| W                     | g                                                    | Water content                               | N/A           | N/A             |
| $D$                   | cm                                                   | Dimensions of compartment                   | [0.5, 5.0]    | None            |
| V                     | cm$^3$                                               | Volume of compartment                       | N/A           | N/A             |
| $\Sigma$F             | g h$^{-1}$                                           | Net water influx                            | N/A           | N/A             |
|                       |                                                      |                                             |               |                 |
| shape                 | /                                                    | Shape of compartment                        | Cylinder      | None            |
| $\phi_D$              | MPa$^{-1}$ h$^{-1}$                                  | Dimensional extensibility                   | 0.02          | DeSchepper[^52] |
| $\epsilon_D$          | MPa                                                  | Dimensional elastic modulus                 | 50            | DeSchepper[^52] |
| $\Gamma$              | MPa                                                  | Yield turgor pressure                       | 0.3           | DeSchepper[^52] |
| $T$                   | K                                                    | Temperature                                 | 298.15        | None            |
| $W_{\text{max}}$      | g                                                    | Water capacity of environmental compartment | 1e6           | None            |
| $W_r$                 | g g$^{-1}$                                           | Relative water content                      | 0.8           | None            |
| $t_{\text{sunrise}}$  | h                                                    | Time of sunrise (hours past midnight)       | 8             | None            |
| $t_{\text{sunset}}$   | h                                                    | Time of sunset (hours past midnight)        | 20            | None            |
| $A_{\text{max}}$      | mol cm$^{-1}$ h$^{-1}$                               | Maximum carbon assimilation rate            | 2e-6          | None            |
| $M_c$                 | h$^{-1}$                                             | Rate of carbon consumption                  | 0.05          | None            |
| $K_s$                 | g h$^{-1}$ MPa$^{-1}$ cm$^{-2}$                      | Specific hydraulic conductivity             | 10            | None            |
| $K$                   | g h$^{-1}$ MPa$^{-1}$                                | Hydraulic conductivity                      | 1000          | DeSchepper[^52] |
| $\eta_{\text{night}}$ | g h$^{-1}$ MPa$^{-1}$ (g h$^{-1}$ MPa$^{-1}$)$^{-1}$ | Relative hydraulic conductivity at night    | 0.1           | Caird[^54]      |
| $P\_h$                | MPa                                                  | Gravitational water potential               | N/A           | None            |
| $h$                   | cm                                                   | Height above ground                         | 0.0           | None            |
