#!/bin/bash
# ===========================================================================
# Master submit - chains the 6-hourly pipeline with Slurm dependencies.
#
#   01 prepare P-files        \
#                       >-- both --> 03 S-files --> 04 LAGRANTO --> 05 decomp --> 06 origin
#   02 T climatology (tbar)   /
#
# Dates, run name, arrival box and levels come from the EVENT block of config.sh.
# The ERA5 period for the P/S files is derived from them:
#   DATA_START = TRAJ_START - 15 days (15-day back-trajectories), DATA_END = TRAJ_END
#
# Usage:
#   ./submit_all.sh                          # everything from config.sh
#   ./submit_all.sh TRAJ_START TRAJ_END      # override the arrival days
#   ./submit_all.sh DATA_START DATA_END TRAJ_START TRAJ_END   # old explicit form
# ===========================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${HERE}"   # stage scripts locate config.sh / *.py via this

case $# in
    0) ;;
    2) export TRAJ_START="$1" TRAJ_END="$2" ;;
    4) DATA_START="$1"; DATA_END="$2"; export TRAJ_START="$3" TRAJ_END="$4" ;;
    *) echo "usage: $0 [TRAJ_START TRAJ_END] | DATA_START DATA_END TRAJ_START TRAJ_END"; exit 1 ;;
esac

source "${HERE}/config.sh"
DATA_START="${DATA_START:-$(date -d "${TRAJ_START} -15 day" +%Y-%m-%d)}"
DATA_END="${DATA_END:-${TRAJ_END}}"
# Slurm logs to /homedata too (the #SBATCH --output default would drop them
# into the code directory you submit from); --chdir keeps any stray relative
# write out of $HOME as well
SB=(sbatch --parsable --partition="${SLURM_PARTITION}"
    --output="${SLURM_LOG_DIR}/%x-%j.out" --chdir="${SLURM_LOG_DIR}")

# ---- pre-flight: fail now rather than after hours in the queue ----
if [[ ! -s "${LAGRANTO_PARAM_DIR}/tracevars.shellscript" ]]; then
    echo "missing ${LAGRANTO_PARAM_DIR}/tracevars.shellscript"; exit 1
fi
if ! grep -Eq '^[[:space:]]*set[[:space:]]+interval[[:space:]]*=[[:space:]]*"?6"?' \
        "${LAGRANTO_ROOT}/goodies/datelist.sh"; then
    echo "set  set interval = \"6\"  in ${LAGRANTO_ROOT}/goodies/datelist.sh first"; exit 1
fi

# LAGRANTO trace on spiritx segfaults with more than ~120 model levels
# (tested: 120 ok, 128 crash) -> refuse anything above 120
if [[ -z "${LEVELS}" ]]; then NLEV=137
elif [[ "${LEVELS}" =~ ^([0-9]+)/([0-9]+)$ ]]; then NLEV=$(( BASH_REMATCH[2] - BASH_REMATCH[1] + 1 ))
else NLEV=$(tr ',' '\n' <<< "${LEVELS}" | grep -c .)
fi
if (( NLEV > 120 )); then
    echo "LEVELS='${LEVELS}' keeps ${NLEV} levels; trace crashes above 120 (use e.g. 18/137)"; exit 1
fi

# climatology base period = data period +/- 1 day (for the +/-6 h neighbours)
BASE_START="$(date -d "${DATA_START} -1 day" +%Y-%m-%d)"
BASE_END="$(date -d "${DATA_END} +1 day" +%Y-%m-%d)"

echo "P/S data period : ${DATA_START} .. ${DATA_END}"
echo "climatology     : ${BASE_START} .. ${BASE_END}"
echo "trajectory runs : ${TRAJ_START} 00 UTC .. ${TRAJ_END} 18 UTC (6-hourly)"
echo "arrival box     : ${ARRIVAL_LON_MIN}..${ARRIVAL_LON_MAX}E x ${ARRIVAL_LAT_MIN}..${ARRIVAL_LAT_MAX}N" \
     "(${ARRIVAL_DLON} x ${ARRIVAL_DLAT} deg), ${ARRIVAL_NLEV} levels ${ARRIVAL_P_BOTTOM}..${ARRIVAL_P_TOP} hPa"
echo "output          : ${NC_BASE}  +  ${BASE_OUT}"

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
