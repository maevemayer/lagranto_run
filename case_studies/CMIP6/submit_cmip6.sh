#!/bin/bash
# ===========================================================================
# CMIP6 LAGRANTO run: prepare_cmip6.py (P/S files) -> 04 LAGRANTO -> 05 decomp
# -> 06 origin, chained with Slurm dependencies. Model, member, dates: config.sh
# (or export CMIP6_MODEL=..., TRAJ_START=..., TRAJ_END=... before calling).
#
# Usage:
#   ./submit_cmip6.sh --test     # one P/S time step + plausibility checks only
#   ./submit_cmip6.sh            # the whole chain
# ===========================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${HERE}"
TEST=0
case "${1:-}" in
    --test) TEST=1 ;;
    "") ;;
    *) echo "usage: $0 [--test]"; exit 1 ;;
esac
source "${HERE}/config.sh"

RUN_TAG="lagranto:${RUN_NAME}"
if squeue -h -u "${USER}" -o "%k" 2>/dev/null | grep -qxF "${RUN_TAG}"; then
    echo "jobs of run ${RUN_NAME} are still queued/running - wait or scancel them first"; exit 1
fi
if ! grep -Eq '^[[:space:]]*set[[:space:]]+interval[[:space:]]*=[[:space:]]*"?6"?' \
        "${LAGRANTO_ROOT}/goodies/datelist.sh"; then
    echo "set  set interval = \"6\"  in ${LAGRANTO_ROOT}/goodies/datelist.sh first"; exit 1
fi
SB=(sbatch --parsable --partition="${SLURM_PARTITION}" --comment="${RUN_TAG}"
    --output="${SLURM_LOG_DIR}/%x-%j.out" --chdir="${SLURM_LOG_DIR}")

echo "model    : ${CMIP6_MODEL} ${CMIP6_EXP} ${CMIP6_MEMBER:-(first member covering the period)} from ${CMIP6_POOL}"
echo "arrivals : ${TRAJ_START} 00 UTC .. ${TRAJ_END} 18 UTC, box ${ARRIVAL_LON_MIN}..${ARRIVAL_LON_MAX}E" \
     "x ${ARRIVAL_LAT_MIN}..${ARRIVAL_LAT_MAX}N, ${ARRIVAL_NLEV} levels ${ARRIVAL_P_BOTTOM}..${ARRIVAL_P_TOP} hPa"
echo "P/S      : ${DATA_DIR} (${GRID_DEG} deg, ${LAT_MIN}-90N)"
echo "results  : ${RESULTS_DIR} (trajectories), ${NC_RESULTS_DIR} (netCDF)"

if (( TEST )); then
    j=$("${SB[@]}" --export=ALL,PREP_ARGS=--test --time=03:00:00 "${HERE}/01_prepare_cmip6.slurm")
    echo "  prepare --test -> ${j}   (log: ${SLURM_LOG_DIR}/prepcmip6-${j}.out)"
    exit 0
fi
j1=$("${SB[@]}" "${HERE}/01_prepare_cmip6.slurm")
echo "  prepare     -> ${j1}"
j4=$("${SB[@]}" --dependency=afterok:${j1} --export=ALL,START="${TRAJ_START}",END="${TRAJ_END}" \
     "${HERE}/04_run_lagranto.slurm")
echo "  04 LAGRANTO -> ${j4}"
j5=$("${SB[@]}" --dependency=afterok:${j4} "${HERE}/05_traj_postprocess.slurm")
echo "  05 decomp   -> ${j5}"
j6=$("${SB[@]}" --dependency=afterok:${j5} "${HERE}/06_origin_contributions.slurm")
echo "  06 origin   -> ${j6}"
echo "Watch with:  squeue -u \$USER"
