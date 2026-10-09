#!/bin/bash
# ===========================================================================
# Master submit - chains the 6-hourly pipeline with Slurm dependencies.
#
#   01 prepare P-files        \
#                       >-- both --> 03 S-files --> 04 LAGRANTO --> 05 decomp --> 06 origin
#   02 T climatology (tbar)   /
#
# DATA_START/DATA_END : period for which P/S files are built. Must cover
#                       [TRAJ_START - 15 days , TRAJ_END] (15-day back traj).
# TRAJ_START/TRAJ_END : the 6-hourly trajectory start times (stage 04).
#
# Usage:
#   ./submit_all.sh DATA_START DATA_END TRAJ_START TRAJ_END
# Example (ERA5T 2026 event, starts 2026-05-10..06-01, 15 d back => data from 2026-04-25):
#   ./submit_all.sh 2026-04-25 2026-06-01 2026-05-10 2026-06-01
# ===========================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${HERE}"   # stage scripts locate config.sh / *.py via this

if [[ $# -ne 4 ]]; then
    echo "usage: $0 DATA_START DATA_END TRAJ_START TRAJ_END"; exit 1
fi
DATA_START="$1"; DATA_END="$2"; TRAJ_START="$3"; TRAJ_END="$4"

source "${HERE}/config.sh"
SB=(sbatch --parsable --partition="${SLURM_PARTITION}")

# ---- pre-flight: fail now rather than after hours in the queue ----
for f in create_startf.param tracevars.shellscript; do
    if [[ ! -s "${LAGRANTO_PARAM_DIR}/${f}" ]]; then
        echo "missing ${LAGRANTO_PARAM_DIR}/${f} (copy it from Levante, see README)"; exit 1
    fi
done
if ! grep -Eq '^[[:space:]]*set[[:space:]]+interval[[:space:]]*=[[:space:]]*"?6"?' \
        "${LAGRANTO_ROOT}/goodies/datelist.sh"; then
    echo "set  set interval = \"6\"  in ${LAGRANTO_ROOT}/goodies/datelist.sh first"; exit 1
fi

# climatology base period = data period +/- 1 day (for the +/-6 h neighbours)
BASE_START="$(date -d "${DATA_START} -1 day" +%Y-%m-%d)"
BASE_END="$(date -d "${DATA_END} +1 day" +%Y-%m-%d)"

echo "P/S data period : ${DATA_START} .. ${DATA_END}"
echo "climatology     : ${BASE_START} .. ${BASE_END}"
echo "trajectory runs : ${TRAJ_START} .. ${TRAJ_END}"

j01=$("${SB[@]}" --export=ALL,START="${DATA_START}",END="${DATA_END}" \
      "${HERE}/01_prepare_era5.slurm")
echo "  01 prepare      -> ${j01}"

j02=$("${SB[@]}" --export=ALL,BASE_START="${BASE_START}",BASE_END="${BASE_END}" \
      "${HERE}/02_t_climatology.slurm")
echo "  02 climatology  -> ${j02}"

j03=$("${SB[@]}" --dependency=afterok:${j01}:${j02} \
      --export=ALL,START="${DATA_START}",END="${DATA_END}" \
      "${HERE}/03_build_s_files.slurm")
echo "  03 S-files       -> ${j03}  (after ${j01} & ${j02})"

j04=$("${SB[@]}" --dependency=afterok:${j03} \
      --export=ALL,START="${TRAJ_START}",END="${TRAJ_END}" \
      "${HERE}/04_run_lagranto.slurm")
echo "  04 LAGRANTO     -> ${j04}  (after ${j03})"

j05=$("${SB[@]}" --dependency=afterok:${j04} \
      "${HERE}/05_traj_postprocess.slurm")
echo "  05 decomposition -> ${j05}  (after ${j04})"

j06=$("${SB[@]}" --dependency=afterok:${j05} \
      "${HERE}/06_origin_contributions.slurm")
echo "  06 origin contrib -> ${j06}  (after ${j05})"

echo "Submitted. Watch with:  squeue -u \$USER"
