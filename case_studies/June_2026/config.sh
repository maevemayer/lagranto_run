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
# Values marked  [CHECK]  are best guesses for spiritx - confirm them with
#   bash tools/check_spiritx_env.sh      (from the repo root)
# before the first submission.
# ===========================================================================

# ---- Python / CDO environment ----------------------------------------------
# Slurm jobs are non-interactive shells, so the conda env is activated
# explicitly. The env needs: python-cdo, xarray, netCDF4, pandas, joblib
# (and a `cdo` binary, either from the same env or from CDO_MODULE).
export CONDA_ENV="${CONDA_ENV:-base}"                         # [CHECK]
export CDO_MODULE="${CDO_MODULE:-}"     # e.g. "cdo/2.2.2" if cdo comes from `module load`  [CHECK]
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

# ---- Slurm ------------------------------------------------------------------
# No account on spiritx. The partition is also hard-coded in the #SBATCH
# headers (they cannot read this file); submit_all.sh passes this value on the
# sbatch command line, which overrides the header.
export SLURM_PARTITION="${SLURM_PARTITION:-zen16}"           # [CHECK] sinfo -s

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
# P/S files are ~0.5 GB each at 137 levels: keep this off /home.
export WORK_ROOT="${WORK_ROOT:-/scratchx/${USER}/lagranto}"    # [CHECK]
export BASE_OUT="${BASE_OUT:-${WORK_ROOT}/2026_May/corrected_6h}"
export SCRATCH_BASE="${SCRATCH_BASE:-${WORK_ROOT}/tmp/2026_MAY_6h}"
export DATA_DIR="${BASE_OUT}/data"        # P- and S-files land here
export TBAR_DIR="${BASE_OUT}/tbar"        # final T climatology (tbar) files
export RESULTS_DIR="${BASE_OUT}/results_6h"
mkdir -p "${DATA_DIR}" "${TBAR_DIR}" "${RESULTS_DIR}" "${SCRATCH_BASE}"

# ---- The 6-hourly cadence ----
export HOURS_CSV="0,6,12,18"              # cdo -selhour argument
export SLOTS="00 06 12 18"               # human-readable slot labels
export STEPS_PER_DAY=4

# ---- Optional vertical restriction (see README) ----
#   empty  -> keep ALL model levels (math identical to notebooks) [DEFAULT]
#   "70/137" -> lowest 68 levels (~250 hPa -> surface); WILL clip 200 hPa
#               trajectory starts and any air above ~250 hPa. Enable only
#               after confirming your trajectories stay below that top.
export LEVELS=""

# ---- ERA5 remap target grid ----
export REMAP_GRID="r720x360"             # 0.5 deg global (unchanged)

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
export LAGRANTO_ROOT="/home/mmayer/code/lagranto.era5.test/lagranto.shellscript"
# tracevars.shellscript used by stage 04 (from Levante)
export LAGRANTO_PARAM_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}/lagranto_params"
# Optional extra library path for the LAGRANTO binaries (netCDF-Fortran), only
# needed if caltra/trace cannot find libnetcdff inside a Slurm job.
export LAGRANTO_LD_LIBRARY_PATH="${LAGRANTO_LD_LIBRARY_PATH:-}"
# IMPORTANT: set  set interval = "6"  in
#   ${LAGRANTO_ROOT}/goodies/datelist.sh
# before running stage 04 (was "3" for the old 3-hourly data).
