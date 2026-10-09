# Corrected 6-hourly LAGRANTO pipeline (IPSL spiritx)

SLURM-submittable version of the notebook workflow, with the **6-hourly
snapshot** correction. All the physics/math is carried over verbatim from the
notebooks; only the items below change.

Originally written for DKRZ Levante; now ported to spiritx (see
[Port from Levante](#port-from-levante-to-spiritx) at the end).

## What changed (and only this)

1. **No 3-hourly means.** `cdo timselmean(3)` is removed from stages 01 and 02.
2. **6-hourly instantaneous snapshots** (00/06/12/18 UTC) are selected with
   `-selhour,0,6,12,18` **before** the spectral transform (`sp2gp`) — fewer
   timesteps through the expensive step, so it's also the main speed-up.
3. **No rename step.** The old `04_02_rename` corrected a +1 h offset that came
   from *centred* 3h-bin means. Snapshots are instantaneous → times are already
   exactly 00/06/12/18 → nothing to rename. (This is also why the separate
   TCLIM +1 h fix is not needed here.)
4. **DTDT stencil = 6 h.** In stage 03 the climatology time derivative uses the
   ±6 h neighbours (`dt_seconds = 2*6*3600`) because the climatology now lives
   on a 6 h grid. Everything else in 03 (THETA, DTDP, PS, assembly) is identical.
5. **Domain cut to the Northern Hemisphere from the start.** Stages 01 and 02
   apply `-sellonlatbox,$DOMAIN_LATBOX` (default `0,360,0,90`) **chained onto the
   remap**, so ERA5 is trimmed the instant it is regridded and *every* downstream
   file (p, dp, P-files, tbar) is NH-only — less scratch I/O and RAM, faster
   `sp2gp`/remap chains. Set `DOMAIN_LATBOX=""` for the old global behaviour, or
   widen the box (e.g. `0,360,-15,90`) for an equator buffer. 15-day back
   trajectories from W Europe stay well inside the NH; any parcel that leaves the
   box is flagged `-999` by LAGRANTO and masked in stage 06.

Nothing else — grid (`r720x360`), 9-year selection, 21-day running mean,
central differences, potential-temperature formula, variable renames — is
untouched.

## Arrival location for land masking

Stage 06 (event **and** analogue) already writes the **exact arrival point**
(`lon0`, `lat0`, `p0`/`arrival_p`) per trajectory alongside TA and every process
contribution, so a land/ocean mask on arrival is a straight
`is_land(lon0, lat0)` filter on `origin_contributions.nc` — for both the ERA5
event and each analogue rank. (The earlier analogue pooled files that lacked
arrival coordinates came from the *old* framework; files produced by this
pipeline carry them.)

## Files

| file | stage |
|---|---|
| `config.sh` | shared paths / cadence / `LEVELS` knob / conda env (source this) |
| `era5_input.py` | reads one field out of the spiritx all-variable daily ERA5 GRIB (used by 01, 02) |
| `lagranto_params/` | `tracevars.shellscript` for stage 04 (`create_startf.param` is generated per run) |
| `01_prepare_era5_6h.py` + `.slurm` | ERA5 → P-files |
| `02_t_climatology_6h.py` + `.slurm` | 9-year T climatology (tbar) |
| `03_build_s_files_6h.py` + `.slurm` | build S-files (hyam/hybm folded in) |
| `04_run_parallel.csh` + `04_run_lagranto_driver.py` + `.slurm` | LAGRANTO |
| `05_traj_postprocess_6h.py` + `.slurm` | trajectory process decomposition |
| `06_origin_contributions_6h.py` + `.slurm` | origin contributions (TA + per-process genesis + exact arrival) |
| `submit_all.sh` | chain 01+02 → 03 → 04 → 05 → 06 with dependencies |

There is **no stage 04**: the old `04_write_coords.sh` (`ncks -A -v hyam,hybm P S`)
is folded into stage 03 — the coefficients are written straight into each S file
in xarray (no second NCO pass, no `module load nco`) — and the `04_02` rename is
obsolete because the 6-hourly snapshots already carry their true 00/06/12/18
timestamps.

**`$HOME` holds code only.** Output is split over two roots (see `config.sh`);
`RUN_NAME` (default `2026_May/corrected_6h`) is the sub-path for one run under both:

```
/scratchx/$USER/lagranto_run/<RUN_NAME>/      NC_ROOT - every netCDF
├── tmp/                             regridding scratch (stages 01, 02)
├── data/                            P- and S-files (~0.5 GB each), LAGRANTO run dirs in data/tmp/
├── tbar/                            T climatology
└── results_6h/                      traj_netcdfs/, origin_contributions.nc

/homedata/$USER/lagranto_run/                 WORK_ROOT - logs + ASCII trajectories
├── slurm_logs/                      %x-%j.out of every job
└── <RUN_NAME>/results_6h/           trajectory_* (LAGRANTO output), logs/
```

**Intermediates** (`tmp/prepare/`, `tmp/t_clim/bar_intermediate/`,
`data/tmp/tmp_*`) are deleted once the stage that used them has succeeded
(`CLEANUP_INTERMEDIATES=1`, the default). Set it to `0` to keep them for
debugging. A rerun of stage 02 skips straight through if all `tbar/` files exist.

**/scratchx may be purged.** Copy `origin_contributions.nc` (and whatever
else you need to keep) to `/homedata` once a run is done; the trajectories on
`/homedata` are enough to rebuild stages 05/06.

Submitting single stages by hand with `sbatch` also logs there (the `#SBATCH
--output` headers point at `/homedata/%u/lagranto_run/slurm_logs/`), but the
directory has to exist: `source config.sh` once creates it.

## Prerequisites (one-time, manual)

- Create the conda env (repo root, once):
  ```bash
  export CONDA_PKGS_DIRS=/homedata/$USER/conda/pkgs
  conda env create -p /homedata/$USER/conda/thesis -f environment.yml
  ```
- Run `bash tools/check_spiritx_env.sh` (repo root) and fix any `[CHECK]` left in
  `config.sh`.
- Set `set interval = "6"` in `$LAGRANTO/goodies/datelist.sh` (was `"3"`).
- `trace` must handle 137 levels: the LAGRANTO source has `parameter (nlevmax=100)`
  in `trace/trace.f90`, which makes `trace` segfault on full ERA5 (caltra is
  fine). Set it to 200 and rebuild `trace` (spiritx: done Oct 2026).

`submit_all.sh` refuses to submit until the datelist interval is 6 and `nlevmax` >= 137.

## Run an event

Everything event-specific is in the **EVENT block at the top of `config.sh`**:

| setting | default (July 2022 heatwave) |
|---|---|
| `RUN_NAME` | `2022_July/corrected_6h` (sub-path under /scratchx and /homedata) |
| `TRAJ_START`, `TRAJ_END` | `2022-07-13`, `2022-07-24`: every 6 h from 13 Jul 00 UTC to 24 Jul 18 UTC |
| `ARRIVAL_LON_MIN/MAX`, `ARRIVAL_LAT_MIN/MAX` | -5.5..5.0 E, 44.5..55.0 N (lon in -180..180, boxes across 0E are fine) |
| `ARRIVAL_DLON`, `ARRIVAL_DLAT` | 0.5, 0.5 deg |
| `ARRIVAL_P_BOTTOM`, `ARRIVAL_P_TOP`, `ARRIVAL_NLEV` | 1000, 200, 17 hPa (every 50 hPa) |
| `CLIM_WINDOW` | `centered`: analogue rule, year +/-4 clamped to 1940-2025 -> 2017..2025 for 2022 |

Then simply:
```bash
./submit_all.sh
```
The ERA5 period is derived: `TRAJ_START - 15 days .. TRAJ_END` (here 2022-06-28 .. 07-24).
`./submit_all.sh 2022-07-15 2022-07-20` overrides just the arrival days.

## Run the 30 best analogues

Each analogue has its own dates/climatology years. Loop over your rank table,
setting the four dates per rank and its own output subdir. Minimal pattern:

```bash
while read RANK DATA_START DATA_END TRAJ_START TRAJ_END; do
  export RUN_NAME="2026_May/corrected_6h/rank${RANK}"
  ./submit_all.sh "$DATA_START" "$DATA_END" "$TRAJ_START" "$TRAJ_END"
done < ranks.txt
```
(`config.sh` derives all subpaths from `RUN_NAME`, so exporting it per rank
keeps each analogue isolated. With the default `CLIM_WINDOW=centered` every
analogue gets its own centred 9-year climatology (year +/-4, clamped to
1940-2025), as in the Levante analogue pipeline.)

## Optional vertical restriction — `LEVELS`

`config.sh` has `LEVELS=""` (all 137 model levels → math identical to the
notebooks). Setting e.g. `LEVELS="70/137"` keeps the lowest 68 levels
(~250 hPa → surface), ~50% less 3D data and a matching regrid speed-up.

**Caveat for this study:** `create_startf` seeds trajectories on
`@profile(ARRIVAL_P_BOTTOM,ARRIVAL_P_TOP,ARRIVAL_NLEV)` (default 1000,200,17) — up to **200 hPa**, above level 70 (~250 hPa) — and
15-day back trajectories often ascend into the UTLS. `70/137` would clip the
200 hPa starts and any air above ~250 hPa. Enable a level cut only after
confirming your trajectories stay below the chosen top (and raise the top, e.g.
~level 55–60 ≈ 150–170 hPa, if not). When `LEVELS` is set it is applied to the
4-D fields only; the p/dp interface-coefficient handling should be validated.

## Rough cost (Levante numbers; spiritx not yet measured)

Anchored to the Levante logs (`pfiles-rank17` 0.64; `tclim-rank20` 1.9 for 6/22 dates
→ ~3.0 fresh full; `lagranto-rank3` 0.22):

- **~1.1–1.4 node-h per analogue** with the hour-selection (this pipeline).
- **30 analogues + event ≈ ~35–45 node-h**; budget ~50 with margin.
- Without hour-selection (regrid all 24 h) it would be ~80–115 node-h.

## Port from Levante to spiritx

What changed for spiritx (the physics is untouched):

| | Levante | spiritx |
|---|---|---|
| ERA5 input | `/pool/data/ERA5/{ET,E5}/{ml,sf}/an/1H/<code>/…_<code>.grb`, one file per variable, spectral T/U/V | `/proju/flexpart/flexpart_in/ERA5/EN-true/YYYY/ERA5ENYYYYMMDD.grb`, all fields in one daily file; picked with `cdo -selname` (`era5_input.py`) |
| grid | native ~0.28° ERA5 → bilinear 0.5° (`r720x360`) | `/proju` ERA5 is 1° (360×181, w + lnsp, 3-hourly); P/S files kept at native 1° (`REMAP_GRID=none`, latitudes flipped south→north), start points still on the 0.5° arrival grid |
| ET/E5 choice | stage 02 used ET for the current year | single pool for all years |
| surface pressure | `sp` (134) | `sp`, or `exp(lnsp)` if `sp` is absent |
| Python/CDO | `module load python3`, spack CDO 2.2.2 | `thesis` conda env from `environment.yml` on `/homedata` (any CDO; also colindex2 + analysis stack) |
| Slurm | `--account=uc1275 --partition=shared`, up to 84 CPUs | no account, `--partition=zen4`, ≤32 CPUs |
| LAGRANTO | `/work/uc1275/…/lagranto.shellscript`, intel/netcdf modules | `/home/mmayer/code/lagranto.era5.test/lagranto.shellscript`, no modules (`LAGRANTO_LD_LIBRARY_PATH` if needed) |
| outputs | `/work/uc1275/…`, `/scratch/u/u301827/…`, Slurm logs in the submit dir | netCDF (incl. P/S) in `/scratchx/$USER/lagranto_run/`, logs + trajectories in `/homedata/$USER/lagranto_run/`, nothing in `$HOME` |
| param files | copied from fixed `/work` paths (incl. a stale `create_startf.param` that create_startf.sh overwrites anyway) | `tracevars.shellscript` versioned in `lagranto_params/` |
| stage 04 logs | discarded (`/dev/null`) | `results_6h/logs/lagranto_<YYYYMMDD_HH>.log`; failed start times are listed and fail the job |

### Open questions about the `/proju` input

Answer these with `tools/check_spiritx_env.sh` before the first run:

1. **Omega.** The pipeline needs `w` (param 135, Pa/s) on model levels. FLEXPART
   input files often carry `etadot` (77) instead. If 135 is missing, the input
   cannot be used as-is.
2. **Domain.** The pipeline expects at least the full NH. If the files are
   regional, the remap leaves missing values outside, and 15-day trajectories
   leave the domain early.
3. **Levels.** 137 levels keeps the math identical. With fewer levels, check that
   the 200 hPa starts and the hybrid coefficients (`hyai/hybi/hyam/hybm` in the
   netCDF that cdo writes) still make sense.
4. **CDO version.** Any version works: `compute_pressure_and_dp` looks up the
   coefficient dims (`nhym`/`nhyi` in cdo 2.2.2) instead of hard-coding them;
   output is bit-identical to the notebook code. Check that the P-files still
   carry `hyam hybm hyai hybi` (LAGRANTO reads them).

Quick one-day test of the ERA5 reader (stage 01) before submitting everything:
```bash
source config.sh
"$PYTHON" 01_prepare_era5_6h.py 2026-05-25 2026-05-25 1
ncdump -h "$DATA_DIR/P20260525_12" | head -40
```
