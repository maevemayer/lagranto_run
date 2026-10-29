# Input and traced variables

Overview of the fields the Lagrangian workflow uses, for ERA5 and CMIP6. LAGRANTO
([Sprenger & Wernli 2015](https://doi.org/10.5194/gmd-8-2569-2015)) reads one *primary* (P) and one
*secondary* (S) netCDF file per 6-hourly time step and computes 15-day backward trajectories; the
temperature anomaly along each trajectory is then decomposed into advection, adiabatic and diabatic
contributions ([Röthlisberger & Papritz 2023](https://doi.org/10.1038/s41561-023-01126-1);
[Hotz et al. 2024](https://doi.org/10.5194/wcd-5-323-2024)).

## Grid, levels and time

| | ERA5 | CMIP6 |
|---|---|---|
| Source | ERA5 model-level data (hourly; DKRZ data pool) | `6hrLev` output (DKRZ data pool / ESPRI mirror) |
| Time steps | 6-hourly instantaneous fields, 00/06/12/18 UTC | 6-hourly instantaneous fields, 00/06/12/18 UTC |
| Horizontal grid | 0.5° × 0.5° (1° in tests), 0–90° N | regular 1° × 1°, 0–90° N, interpolated from the native grid after the vertical velocity is computed |
| Vertical levels | native hybrid σ–p model levels 18–137 (120 levels, top ≈ 1 hPa) | all native hybrid σ–p model levels (e.g. 80 for MRI-ESM2-0) |

## Primary files (`PYYYYMMDD_HH`)

| Variable | Description | Units | ERA5 field | CMIP6 field |
|---|---|---|---|---|
| `T` | temperature | K | `t` | `ta` |
| `U`, `V` | horizontal wind components | m s⁻¹ | `u`, `v` | `ua`, `va` |
| `OMEGA` | vertical velocity in pressure coordinates | Pa s⁻¹ | `w` | **diagnosed** from `ua`, `va`, `ps` (mass continuity in hybrid coordinates, [Simmons & Burridge 1981](https://doi.org/10.1175/1520-0493(1981)109%3C0758:AEAAMC%3E2.0.CO;2)); validated against each model's own ω (r = 0.92–1.00 in the troposphere) |
| `PS` | surface pressure | hPa | `sp` (or `exp(lnsp)`) | `ps` |
| `hyam`, `hybm`, `hyai`, `hybi` | hybrid coefficients of the levels in the file; pressure p = `hyam` + `hybm` · `PS` | Pa, 1 | from the files | from the files |
| `p`, `dp` | full-level pressure and layer thickness (ERA5 only; diagnostic, not read by LAGRANTO) | Pa | computed | – |

## Secondary files (`SYYYYMMDD_HH`)

| Variable | Description | Units |
|---|---|---|
| `THETA` | potential temperature, θ = T (1000 hPa / p)^0.286 | K |
| `TCLIM` | climatological temperature T̄: 21-day running mean at the same 6-hourly time slot, averaged over a centred 9-year window around the event year (the dataset's own climatology: ERA5 for ERA5, each model for itself) | K |
| `DTDT` | its time tendency ∂T̄/∂t, centred difference over ±6 h | K s⁻¹ |
| `PS`, `hyam`, `hybm` | as in the primary files | |
| `Q` | specific humidity, for the moisture source and sink diagnostics (*planned*; ERA5 `q`, CMIP6 `hus`) | kg kg⁻¹ |

## Variables traced along each trajectory

`T`, `OMEGA`, `THETA`, `TCLIM`, `TCLIM` 10 hPa below and above the parcel (for the vertical gradient
∂T̄/∂p along the trajectory, as in Hotz et al. 2024) and `DTDT`
([`lagranto_params/tracevars.shellscript`](case_studies/June_2026/lagranto_params/tracevars.shellscript)).

## Trajectory set-up (case studies)

- Start points: 0.5° grid over the heatwave region (here 5.5° W–5° E, 44.5–55° N), every 50 hPa from
  1000 to 200 hPa (17 levels), i.e. about 8 000 trajectories per start time
- Start times: every 6 h over the heatwave period (12 days → about 400 000 trajectories per event)
- Length: 15 days backward

## Temperature anomaly decomposition

With the anomaly T′ = T − T̄ and ω = dp/dt along the trajectory, each 6-hour step is split into

| Term | Expression |
|---|---|
| advection | − dT̄/dt + ∂T̄/∂t + ω ∂T̄/∂p |
| adiabatic | (κ T / p − ∂T̄/∂p) ω |
| diabatic | (p / p₀)^κ dθ/dt |
| seasonal | − ∂T̄/∂t |

(κ = R/c_p, p₀ = 1000 hPa; d/dt along the trajectory, ∂/∂t at a fixed point.) The terms are integrated
from the genesis of the anomaly to the arrival point. Code: stages 01–06 in
[`case_studies/June_2026`](case_studies/June_2026) (ERA5) and
[`case_studies/CMIP6`](case_studies/CMIP6) (CMIP6 input files, incl. the vertical velocity in
[`tools/omega`](tools/omega)).
