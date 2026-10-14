#!/bin/bash
# ===========================================================================
# Shared configuration for the corrected 6-hourly LAGRANTO pipeline (IPSL spiritx).
# Source this from every stage:  source config.sh
#
# CORRECTION vs the old notebooks:
#   * NO 3-hourly means (cdo timselmean removed everywhere).
#   * The 6-hourly INSTANTANEOUS snapshots (00/06/12/18 UTC) are selected
#     directly from the hourly ERA5 files, BEFORE any regridding.
#   * Because snapshots are instantaneous (not centred 3h-bin means), the
#     old +1 h timestamp offset is gone -> the 04_02 rename step is dropped.
#   * The DTDT climatology time-stencil is 6 h instead of 3 h (forced by the
#     new cadence). All other math is unchanged from the notebooks.
#
# New user? Follow "Setting up" in the top-level README, then check with
#   bash tools/check_spiritx_env.sh      (from the repo root)
# ===========================================================================

# ===================== EVENT - edit this block per run =======================
# Every value can also be overridden from the shell (export X=... before
# ./submit_all.sh), e.g. in a loop over analogues.

# Run name: sub-path for this run under /scratchx and /homedata
export RUN_NAME="${RUN_NAME:-2022_July/corrected_6h}"

# Arrival (trajectory start) days, inclusive: every 00/06/12/18 UTC from
# TRAJ_START 00 UTC to TRAJ_END 18 UTC. ERA5 is prepared automatically from
# TRAJ_START - 15 days (15-day back-trajectories) to TRAJ_END.
export TRAJ_START="${TRAJ_START:-2022-07-13}"
export TRAJ_END="${TRAJ_END:-2022-07-24}"

# Arrival region: a regular lon/lat grid of start points (degrees; longitudes
# in -180..180 so boxes across 0E work, e.g. -5.5 .. 5.0).
export ARRIVAL_LON_MIN="${ARRIVAL_LON_MIN:--5.5}"
export ARRIVAL_LON_MAX="${ARRIVAL_LON_MAX:-5.0}"
export ARRIVAL_LAT_MIN="${ARRIVAL_LAT_MIN:-44.5}"
export ARRIVAL_LAT_MAX="${ARRIVAL_LAT_MAX:-55.0}"
export ARRIVAL_DLON="${ARRIVAL_DLON:-0.5}"
export ARRIVAL_DLAT="${ARRIVAL_DLAT:-0.5}"

# Arrival pressure levels: ARRIVAL_NLEV levels evenly spaced from
# ARRIVAL_P_BOTTOM up to ARRIVAL_P_TOP (hPa) at every start point
# (LAGRANTO @profile(bottom,top,n)). 1000,200,17 -> every 50 hPa.
export ARRIVAL_P_BOTTOM="${ARRIVAL_P_BOTTOM:-1000}"
export ARRIVAL_P_TOP="${ARRIVAL_P_TOP:-200}"
export ARRIVAL_NLEV="${ARRIVAL_NLEV:-17}"

# Climatology (TCLIM, stage 02) - the paper definition (Hotz et al. 2024): the
# 21-day running-mean T at each 6-hourly slot, averaged over a centred 9-year
# window (event year +/- 4), clamped to 1940..2025 so it never extends past
# 2025: 2021-2026 all get 2017-2025, 2022 -> 2017-2025, 1990 -> 1986-1994.
# CLIM_YEARS="1991-2020" or "2014 2016 2018" overrides it. Years after
# CLIM_LAST_YEAR are refused; stage 02 stops if a year has no ERA5 input
# (/proju starts in 1979). Check finished files: tools/check_climatology_years.py
export CLIM_WINDOW="${CLIM_WINDOW:-centered}"   # only "centered" (paper rule) is accepted
export CLIM_NYEARS="${CLIM_NYEARS:-9}"
export CLIM_FIRST_YEAR="${CLIM_FIRST_YEAR:-1940}"
export CLIM_LAST_YEAR="${CLIM_LAST_YEAR:-2025}"
export CLIM_YEARS="${CLIM_YEARS:-}"
# =============================================================================

# ---- Python / CDO environment ----------------------------------------------
# Slurm jobs are non-interactive shells, so the conda env is activated
# explicitly. The env needs: python-cdo, xarray, netCDF4, pandas, joblib
# (and a `cdo` binary, either from the same env or from CDO_MODULE).
export CONDA_ENV="${CONDA_ENV:-/homedata/${USER}/conda/thesis}"   # from environment.yml (repo root)
export CDO_MODULE="${CDO_MODULE:-}"     # cdo comes from the conda env; /usr/bin/cdo is 1.9.9 (too old)
export CDO_BIN_DIR="${CDO_BIN_DIR:-}"   # optional: directory prepended to PATH for cdo

_cfg_flags="$-"; set +u          # module/conda scripts trip over `set -u`
if [[ -n "${CDO_MODULE}" ]] && command -v module >/dev/null 2>&1; then
    module load "${CDO_MODULE}"
fi
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
fi
if [[ "${_cfg_flags}" == *u* ]]; then set -u; fi
if [[ -n "${CDO_BIN_DIR}" ]]; then export PATH="${CDO_BIN_DIR}:${PATH}"; fi
PYTHON="$(command -v python3 || command -v python)"
export PYTHON
export PYTHONDONTWRITEBYTECODE=1        # no __pycache__ in the code tree

# ---- Slurm ------------------------------------------------------------------
# No account on spiritx. The partition is also hard-coded in the #SBATCH
# headers (they cannot read this file); submit_all.sh passes this value on the
# sbatch command line, which overrides the header.
export SLURM_PARTITION="${SLURM_PARTITION:-zen4}"   # 64 cores / 254 GB nodes (4 GB/core)

# ---- ERA5 input (FLEXPART-input pool on /proju) -----------------------------
# One GRIB file per day holding ALL fields (model-level T/U/V/W + surface):
#   ${ERA5_DIR}/YYYY/ERA5ENYYYYMMDD.grb
export ERA5_DIR="/proju/flexpart/flexpart_in/ERA5/EN-true"
export ERA5_FILE_PATTERN="{year}/ERA5EN{year}{month:02d}{day:02d}.grb"
# GRIB shortNames used for each pipeline field (cdo -selname). Surface pressure
# falls back to exp(lnsp) automatically if `sp` is not in the file.
export ERA5_NAME_T="t"
export ERA5_NAME_U="u"
export ERA5_NAME_V="v"
export ERA5_NAME_W="w"            # omega (param 135). See README if only etadot (77) exists.
export ERA5_NAME_SP="sp"

# ---- Output tree ------------------------------------------------------------
# $HOME holds code only. The pipeline writes to two roots:
#   NC_ROOT   (/scratchx) : every netCDF - P/S files (~0.5 GB each), tbar,
#                           regridding temp files, traj_netcdfs/,
#                           origin_contributions.nc, LAGRANTO run dirs
#   WORK_ROOT (/homedata) : Slurm logs, ASCII trajectory_* files, stage 04 logs
# RUN_NAME (EVENT block above) is the sub-path for one run under both roots.
export NC_ROOT="${NC_ROOT:-/scratchx/${USER}/lagranto_run}"
export WORK_ROOT="${WORK_ROOT:-/homedata/${USER}/lagranto_run}"
export NC_BASE="${NC_ROOT}/${RUN_NAME}"
export BASE_OUT="${WORK_ROOT}/${RUN_NAME}"
export SCRATCH_BASE="${NC_BASE}/tmp"      # regridding scratch (stages 01, 02)
export DATA_DIR="${NC_BASE}/data"         # P- and S-files land here
export TBAR_DIR="${NC_BASE}/tbar"         # final T climatology (tbar) files
export NC_RESULTS_DIR="${NC_BASE}/results_6h"   # traj_netcdfs/, origin_contributions.nc
export RESULTS_DIR="${BASE_OUT}/results_6h"     # trajectory_* (ASCII), logs/
export SLURM_LOG_DIR="${WORK_ROOT}/slurm_logs"  # same dir as the #SBATCH --output headers
# Delete intermediates once the stage that used them has succeeded:
#   01 -> tmp/prepare/ day files, 02 -> tmp/t_clim/bar_intermediate/ running
#   means, 04 -> data/tmp/tmp_<YYYYMMDD_HH>/ run folders. Each stage checks its
#   own outputs first; anything from a failed day/start time is kept.
#   0 keeps everything (debugging; a rerun then skips work already done).
export CLEANUP_INTERMEDIATES="${CLEANUP_INTERMEDIATES:-1}"
mkdir -p "${DATA_DIR}" "${TBAR_DIR}" "${NC_RESULTS_DIR}" "${RESULTS_DIR}" "${SCRATCH_BASE}" "${SLURM_LOG_DIR}"

# ---- The 6-hourly cadence ----
export HOURS_CSV="0,6,12,18"              # cdo -selhour argument
export SLOTS="00 06 12 18"               # human-readable slot labels
export STEPS_PER_DAY=4

# ---- Model levels kept (cdo -sellevidx, see README) ----
#   "18/137" -> 120 levels, model top ~1 hPa [DEFAULT on spiritx]. The spiritx
#               LAGRANTO trace segfaults with more than ~120 levels (120 ok,
#               128 crash); levels 1-17 lie far above any trajectory, so results
#               below are identical to a 137-level run, with ~12% less data.
#   ""       -> all 137 levels (crashes trace on spiritx).
#   "70/137" -> lowest 68 levels (~250 hPa): would clip the 200 hPa starts.
export LEVELS="${LEVELS:-18/137}"

# ---- Grid of the P/S files ----
#   "none"     -> keep the input grid: the /proju ERA5 is 1 deg (360x181), so
#                 P/S files are 1 deg; latitudes are flipped to south->north
#                 (no interpolation). [DEFAULT on spiritx]
#   "r720x360" -> bilinear to 0.5 deg as on Levante (4x the size; from 1 deg
#                 input this adds no real resolution).
export REMAP_GRID="${REMAP_GRID:-none}"

# ---- Domain cut (applied right after remap, in stages 01 & 02) ----
#   ERA5 is cut to this lon,lat box the moment it is regridded, so every
#   downstream file (p, dp, P-files, tbar) is domain-only -> less I/O and RAM.
#   Format: lon1,lon2,lat1,lat2 (CDO -sellonlatbox).
#     "0,360,0,90"   -> whole Northern Hemisphere [DEFAULT]
#     "0,360,-15,90" -> NH plus an equator buffer
#     empty          -> keep the global grid (no cut).
#   If the /proju input is itself regional, everything outside it is missing
#   after the remap; trajectories leaving it are flagged -999 (see README).
export DOMAIN_LATBOX="0,360,0,90"

# ---- LAGRANTO ---------------------------------------------------------------
# your own LAGRANTO build (see README "Setting up"); override with export LAGRANTO_ROOT=...
export LAGRANTO_ROOT="${LAGRANTO_ROOT:-${HOME}/code/lagranto.era5.test/lagranto.shellscript}"
# tracevars.shellscript used by stage 04 (from Levante)
export LAGRANTO_PARAM_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}/lagranto_params"
# Optional extra library path for the LAGRANTO binaries (netCDF-Fortran), only
# needed if caltra/trace cannot find libnetcdff inside a Slurm job.
export LAGRANTO_LD_LIBRARY_PATH="${LAGRANTO_LD_LIBRARY_PATH:-}"
# IMPORTANT: set  set interval = "6"  in
#   ${LAGRANTO_ROOT}/goodies/datelist.sh
# before running stage 04 (was "3" for the old 3-hourly data).
