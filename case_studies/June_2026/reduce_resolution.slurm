#!/bin/bash
#SBATCH --job-name=reduce6h
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/homedata/%u/lagranto_run/slurm_logs/%x-%j.out
# ---------------------------------------------------------------------------
# Reduced-resolution P/S files for a resolution experiment (reduce_resolution.py).
# Normally submitted by submit_resolution_experiment.sh, which sets
#   SRC_DATA_DIR, DST_DATA_DIR, NLEV and (optional) RES.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env (cdo, python), sets PYTHON

: "${SRC_DATA_DIR:?}" "${DST_DATA_DIR:?}" "${NLEV:?}"
"${PYTHON}" "${HERE}/reduce_resolution.py" "${SRC_DATA_DIR}" "${DST_DATA_DIR}" --nlev "${NLEV}" \
    ${RES:+--res "${RES}"} --jobs "${SLURM_CPUS_PER_TASK:-16}"
