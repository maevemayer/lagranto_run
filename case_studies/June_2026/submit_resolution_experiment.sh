#!/bin/bash
# ===========================================================================
# Resolution experiment: reduced copy of an existing run's P/S files, then
# LAGRANTO (04), decomposition (05) and origin contributions (06) as a NEW,
# separate run. The source run is not touched.
#
#   NLEV : 60 or 30 model levels below ~100 hPa (reduce_resolution.py)
#   RES  : horizontal grid in degrees (1.5, 2, ...); omit to keep the grid
#
# The new run is  <dir of BASE_RUN>/L<NLEV>_<RES>deg  (RES 1.5 -> 1p5; no RES
# -> 1deg, the native spiritx grid), e.g. 2022_July/L30_2deg, with its own
# data/ on /scratchx and results_6h/ on /homedata. Arrival box, levels and
# times are those of config.sh (same start points as the source run).
#
# Usage:
#   ./submit_resolution_experiment.sh BASE_RUN NLEV [RES]
#   ./submit_resolution_experiment.sh 2022_July/corrected_6h 60
#   ./submit_resolution_experiment.sh 2022_July/corrected_6h 30 2
# ===========================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${HERE}"

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "usage: $0 BASE_RUN NLEV [RES]"; exit 1
fi
BASE_RUN="$1"; NLEV="$2"; RES="${3:-}"
case "${NLEV}" in 60|30) ;; *) echo "NLEV must be 60 or 30"; exit 1 ;; esac
TAG="${RES:-1}"; TAG="${TAG/./p}"
export RUN_NAME="$(dirname "${BASE_RUN}")/L${NLEV}_${TAG}deg"
source "${HERE}/config.sh"            # paths of the NEW run (RUN_NAME set above)

export SRC_DATA_DIR="${NC_ROOT}/${BASE_RUN}/data"
export DST_DATA_DIR="${DATA_DIR}"
export NLEV RES
if ! ls "${SRC_DATA_DIR}"/P[0-9]* >/dev/null 2>&1; then
    echo "no P-files in ${SRC_DATA_DIR}"; exit 1
fi
RUN_TAG="lagranto:${RUN_NAME}"
if squeue -h -u "${USER}" -o "%k" 2>/dev/null | grep -qxF "${RUN_TAG}"; then
    echo "jobs of run ${RUN_NAME} are still queued/running - wait or scancel them first"; exit 1
fi
SB=(sbatch --parsable --partition="${SLURM_PARTITION}" --comment="${RUN_TAG}"
    --output="${SLURM_LOG_DIR}/%x-%j.out" --chdir="${SLURM_LOG_DIR}")

echo "experiment ${RUN_NAME}: ${NLEV} levels, grid ${RES:-unchanged (1 deg)}"
echo "  from ${SRC_DATA_DIR}"
echo "  P/S  ${DST_DATA_DIR}"
echo "  results ${RESULTS_DIR} (trajectories), ${NC_RESULTS_DIR} (netCDF)"
echo "  arrivals ${TRAJ_START} .. ${TRAJ_END}"

j0=$("${SB[@]}" "${HERE}/reduce_resolution.slurm")
echo "  reduce    -> ${j0}"
j4=$("${SB[@]}" --dependency=afterok:${j0} --export=ALL,START="${TRAJ_START}",END="${TRAJ_END}" \
     "${HERE}/04_run_lagranto.slurm")
echo "  04 LAGRANTO -> ${j4}"
j5=$("${SB[@]}" --dependency=afterok:${j4} "${HERE}/05_traj_postprocess.slurm")
echo "  05 decomp   -> ${j5}"
j6=$("${SB[@]}" --dependency=afterok:${j5} "${HERE}/06_origin_contributions.slurm")
echo "  06 origin   -> ${j6}"
echo "Watch with:  squeue -u \$USER"
