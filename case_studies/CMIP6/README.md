# CMIP6 case: LAGRANTO on CMIP6 6hrLev output (spiritx)

Same trajectories and decomposition as `June_2026`, but on model output from the ESPRI
CMIP6 mirror (`/bdd/CMIP6`):

| step | script | |
|---|---|---|
| P/S files | `prepare_cmip6.py` (`01_prepare_cmip6.slurm`) | replaces ERA5 stages 01–03 |
| 04–06 | links to `../June_2026/` | unchanged |

`prepare_cmip6.py`
- reads 6hrLev `ua, va, ps, ta` of one member;
- diagnoses ω from the continuity equation on the native grid with the validated functions of
  `tools/omega/omega_check.py` (spherical harmonics on Gaussian grids, Simmons–Burridge full levels);
- builds TCLIM from the model's own `ta` (21-day running mean per 6-h slot, centred 9-year window
  inside 1950–2014) and DTDT = ±6 h centred difference;
- interpolates to a regular 1° grid 0–90N, writes P (T, U, V, OMEGA, PS) and S (THETA, TCLIM, DTDT, PS)
  files with lev = 1..N and the coefficients of those N levels.

Not yet supported: regridded models (CNRM family, EC-Earth3 → need `omega_check_cnrm.py --unregrid`
first) and noleap/360-day calendars.

## Run

```bash
cd ~/code/lagranto_run/case_studies/CMIP6
./submit_cmip6.sh --test      # one time step + checks (ω std, T − TCLIM, PS range) in the prepcmip6 log
./submit_cmip6.sh             # prepare -> 04 -> 05 -> 06
./submit_cmip6.sh --from-04   # P/S files already there: 04 -> 05 -> 06
```

Model, member, dates, grid: `config.sh` or exported variables, e.g.
`TRAJ_START=2005-08-01 TRAJ_END=2005-08-03 ./submit_cmip6.sh`.
Output: `/scratchx/$USER/lagranto_run/CMIP6/<model>_<exp>/<TRAJ_START>/` and the same path on `/homedata`.
