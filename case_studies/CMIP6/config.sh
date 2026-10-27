#!/bin/bash
# ===========================================================================
# CMIP6 case: the 6-hourly LAGRANTO pipeline on CMIP6 6hrLev model output.
# Stages 01-03 (ERA5) are replaced by prepare_cmip6.py; stages 04-06 are the
# unchanged June_2026 scripts (linked into this folder). Everything not set
# here comes from ../June_2026/config.sh (env, paths, Slurm, LAGRANTO).
# ===========================================================================

# ===================== MODEL + EVENT - edit per run ==========================
export CMIP6_MODEL="${CMIP6_MODEL:-MRI-ESM2-0}"   # native Gaussian grid; regridded models are refused
export CMIP6_EXP="${CMIP6_EXP:-historical}"
export CMIP6_MEMBER="${CMIP6_MEMBER:-}"            # empty: first member with ua/va/ps/ta for the period
export CMIP6_POOL="${CMIP6_POOL:-/bdd/CMIP6}"      # ESPRI mirror (DKRZ pool layout)

# Arrival days (prototype event: replace by a TXx date of the model later)
export TRAJ_START="${TRAJ_START:-2010-07-15}"
export TRAJ_END="${TRAJ_END:-2010-07-17}"

# Model climatology: centred 9-year window, kept inside the run's years
export CLIM_NYEARS="${CLIM_NYEARS:-9}"
export CLIM_FIRST_YEAR="${CLIM_FIRST_YEAR:-1950}"
export CLIM_LAST_YEAR="${CLIM_LAST_YEAR:-2014}"

# Grid of the P/S files (linear from the native grid), Northern Hemisphere
export GRID_DEG="${GRID_DEG:-1.0}"
export LAT_MIN="${LAT_MIN:-0}"

export RUN_NAME="${RUN_NAME:-CMIP6/${CMIP6_MODEL}_${CMIP6_EXP}${CMIP6_MEMBER:+_${CMIP6_MEMBER}}/${TRAJ_START}}"
# =============================================================================

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../June_2026/config.sh"
