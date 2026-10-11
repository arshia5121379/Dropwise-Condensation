# Dropwise Condensation Heat Transfer – Piecewise Model (MATLAB)

`Dropwise_Condensation_HT_Coefficient_PieceWise.m`

## Short description

A MATLAB script that calculates the **heat transfer coefficient** and the **heat transfer rate** for
**dropwise condensation of water vapour from a humid-air stream** flowing along a cooled pipe/channel wall.
The wall is divided into short segments. In each segment the script solves the vapour–air interface
temperature iteratively, combines the droplet-side, gas-side (vapour diffusion with non-condensable air) and
sensible-convection resistances into a total coefficient, computes the heat removed, and updates the bulk gas
temperature before moving to the next segment. It reports the heat transfer per segment, the total heat
transfer, the outlet temperature, and plots heat transfer and bulk temperature along the length.

---

## 1. Requirements

- MATLAB on **Windows** with **NIST REFPROP** installed (the script calls `refpropm` for all fluid properties).
- The REFPROP folder is added at the top of the script. Change this line if yours is elsewhere:
  ```matlab
  addpath('C:\Program Files (x86)\REFPROP');
  ```
- Fluids used: `WATER` and `AIR.MIX`.

## 2. How to run

1. Open the script and edit the **INPUTS** section (Section 3).
2. Run it. The console prints, for every segment, the inlet and outlet bulk temperature, the heat flux and the
   heat transfer rate, followed by the final totals.
3. A figure with two plots opens: heat transfer vs. position, and bulk temperature vs. position.

## 3. Inputs

All temperatures are converted to kelvin inside the script.

| Variable | Default | Meaning |
|---|---|---|
| `T_W` | 5 °C | Wall temperature [K] |
| `T_B` | 25 °C | Bulk (inlet) temperature of the humid-air stream [K] |
| `T_exit` | `T_B - 2` | Target exit temperature [K] (defined but not used in the calculation) |
| `RH` | 1.0 | Relative humidity of the inflow [–] (1 = saturated) |
| `D_h` | 0.0158 | Hydraulic diameter [m] |
| `P_amb` | 1e5 | Total pressure [Pa] |
| `theta` | 150° | Contact angle [rad] (defined but not used; the advancing angle drives the model) |
| `theta_A` | 155° | Advancing contact angle [rad] (`theta_max`) |
| `Diff_droplet` | 2.9482e-5 | Diffusion coefficient used in the Schmidt number [m²/s] |
| `g` | 9.81 | Gravitational acceleration [m/s²] |
| `r_d_mean` | 0.5e-3 | Mean droplet radius [m] (starting value of the maximum droplet radius) |
| `x_total` | 0.1 | Total condensation length [m] |
| `dx` | 0.005 | Segment length [m]. Number of segments = `round(x_total/dx)` |
| `Alpha` | 90° | Inclination angle of the condensing surface [rad] (90° = vertical) |
| `total_HT_area` | 0.062368 | Total heat transfer area [m²]; each segment gets `dA = total_HT_area/n_segments` |
| `A1` | π·(0.035)² | Upstream flow area [m²] (used only to get the channel velocity) |
| `A2` | 0.0024261 | Channel cross-sectional area [m²] |
| `U1` | 2.5 | Upstream velocity [m/s] |
| `U_in` | `A1*U1/A2` | Channel inlet velocity from continuity [m/s] |

Other numbers fixed inside the loop: initial interface-temperature guess `T_ik = 280 K`, convergence tolerance
`eps = 1e-5`, maximum 100 iterations, surface-enhancement constants `m = 1.451`, `n = 0.933` (from the
referenced paper's table 4).

## 4. What the script calculates

### 4.1 Bulk properties (once, at `T_B`)
Latent heat, thermal diffusivity, specific heats, conductivities, surface tension, densities, kinematic and
dynamic viscosities and molar masses of water, vapour and air. From `RH` and the saturation pressure it gets the
vapour and air partial pressures, the **mole fractions** (`Y_AIR`, `Y_H2O`), the mixture molar mass and the
**mass fractions** (`X_AIR`, `X_H2O`). Mixture density and viscosity (`Rho_mix`, `Mu_mix`) are mole-fraction
weighted.

### 4.2 Segment loop (`i = 1 … n_segments`)
For each segment, at position `x_current = i*dx`, the interface temperature `T_i` is found by fixed-point
iteration. Each iteration does:

1. **Interface composition** – saturation pressure of vapour at `T_i` gives the interface mole fraction of
   vapour and air, then the mass fraction of air `X_AIR_i` (saturation-line assumption).
2. **Droplet-side coefficient `h_d`** – from the Reynolds number and the mass-fraction ratio
   `B_i = ln(X_AIR_i / X_AIR)`.
3. **Maximum droplet radius `R_max`** – droplet shape from the advancing/receding contact angles
   (`theta_max`, `theta_min` corrected with the Bond number), drag coefficient `C_drag`, and a balance of
   gravity, drag and surface-tension retention force. It is solved as the smallest root magnitude of the
   polynomial `a²R⁸ + b²R⁴ − c² = 0`.
4. **Surface enhancement factor** `C_Film_Drop = m*R_max + n`, which scales the film-condensation Nusselt number
   to account for droplets.
5. **Dimensionless numbers** – friction factor `ksi`, `Sc`, `Pr`, Jakob number `Ja`, and entrance-corrected
   `Sh0`, `Nu0` (Gnielinski form); Lewis number `Le = (Nu0/Sh0)^(1/3)`.
6. **Correction coefficients** for heat and mass transfer: suction (`SIGMA_sq`, `SIGMA_sm`) and fog formation
   (`SIGMA_fq`, `SIGMA_fm`), using the slope of the saturation line `dF_dT`.
7. **Gas-side coefficient `hg`** from the modified Nusselt number `Nu_g` (heat transfer plus the latent term
   `ln(X_AIR_i/X_AIR)`), mixture conductivity and `D_h`.
8. **Sensible convection `h_convection`** from `Nu = 0.0265·Re^0.8·Pr^0.3` (Incropera, Eq. 8-60 as referenced).
9. **Coefficient combination**
   - `h_t = (1/hg + 1/h_d)^-1`  (condensation, series resistances)
   - `h_Total = h_t + h_convection`
10. **Heat flux** `q_Total = h_Total·(T_B_current − T_W)`, `q_g = q_Total − q_convection`, and the updated
    interface temperature `T_ik = T_B_current − q_g/hg`.

The loop stops when `|T_ik − T_i_new| ≤ eps` or after 100 iterations.

### 4.3 Segment energy balance
- Segment heat transfer: `q_i(i) = q_Total · dA` [W]
- Mass flow: `m_dot_a = Rho_mix · U_in · A2`
- Outlet temperature: `T_out = T_B_current − q_Total·dA / (m_dot_a · cp_mix · 1000)`
- The outlet temperature becomes the next segment's bulk temperature.

## 5. Outputs

| Item | Description |
|---|---|
| Console, per segment | Position, inlet/outlet bulk temperature, heat flux `q_Total` [W/m²], heat transfer `Q_i` [W] |
| `total_heat_transfer` | Sum of all segment heat transfers [W] |
| Final outlet temperature, temperature drop | From `T_out_i(end)` and `T_B` |
| `q_i`, `T_out_i` | Arrays of segment heat transfer [W] and outlet temperature [K] |
| Figure | Subplot 1: heat transfer vs. position. Subplot 2: bulk temperature vs. position |

The heat transfer **coefficients** of the last segment remain in the workspace after the run (`h_d`, `hg`,
`h_convection`, `h_t`, `h_Total`). They are not stored per segment; add arrays (like `q_i`) if you need the
coefficient along the length.

## 6. Assumptions

- Steady flow; wall at a uniform, constant temperature.
- Bulk properties, humidity and mixture composition are evaluated once at the inlet bulk temperature and kept
  constant. Only the bulk temperature is updated segment by segment (the water removed is not subtracted from
  the flow).
- The vapour–air interface is saturated at the interface temperature.
- Dropwise condensation is represented as a surface enhancement factor on a film-condensation Nusselt number.
- Spherical-cap droplets with one mean droplet radius as the starting value.

## 7. Items to check in the code

These are points I noticed while reading the script. They may affect the numbers, so please confirm them:

1. **`A1` and `A2` are overwritten inside the loop.** They first hold the flow areas, but the droplet-shape step
   reassigns them to dimensionless coefficients (`A1 = 2*(cos(theta_max)-cos(theta_min))/pi^3`, `A2 = -3*…`). The
   mass flow `m_dot_a = Rho_mix*U_in*A2` therefore uses that coefficient and not the channel area, which changes
   the temperature drop. Rename the droplet-shape coefficients (for example `Ad1`, `Ad2`) and keep `A2` as the
   channel area.
2. **`sin(beta)` in the drag-term `b`.** `beta = 1 + 0.096*Bo` is a Bond-number factor, not an angle. The
   expression may have been meant to use a contact angle such as `theta_min`.
3. **Units of `h_d`.** It is computed as `0.347·Re_B / (1 + 1.75·exp(−330.5·B_i))`, a dimensionless group, but it
   is combined with `hg` [W/m²K] as a coefficient. Check the reference for the intended dimensional form.
4. **`X_vapor_i` uses the bulk `M_total`** rather than the interface mixture molar mass.
5. **No convergence warning.** If the interface loop reaches 100 iterations, the script continues silently.
6. **Unused variables:** `T_exit`, `theta`, `q_t`, `q_d`, `Alpha_WATER`, `Nu_WATER`, `Nu_mix`, `Alpha_mix`, the
   `h_fg_i` values and `Re_D`'s other uses. Safe to remove or put to use.
7. The `Debug` comments on the `Sh0` and `Nu0` lines and the commented-out fixed values (47.8, 57.5) suggest
   these were being tested. Make sure the active expressions are the intended ones.

## 8. References named in the code

The code comments refer to equation and table numbers (Eq. 2, 8-60, 18, 24, 37, tables 2 and 4) from the work on
turbulent vapour condensation with non-condensable gases and fog formation ("Effect of Fog Formation on Turbulent
Vapor Condensation With Noncondensable Gases") and from Incropera's heat transfer text. Check the exact
equations against those sources.
