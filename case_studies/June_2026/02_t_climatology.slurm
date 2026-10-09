#!/bin/bash
#SBATCH --job-name=tclim6h
#SBATCH --partition=zen16           # [CHECK] spiritx partition (sinfo -s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
# ---------------------------------------------------------------------------
# Stage 02: 9-year T climatology (tbar), 6-hourly.
# BASE_START/END should extend +/-1 day beyond the S-file period (stage 03
# needs the +/-6 h climatology neighbours across day boundaries).
# Submit:  sbatch --export=ALL,BASE_START=2026-05-14,BASE_END=2026-06-02 02_t_climatology.slurm
# ---------------------------------------------------------------------------
set -euo pipefail
# Slurm runs a spool copy of this script, so BASH_SOURCE is useless here:
# use PIPELINE_DIR (exported by submit_all.sh), else the sbatch submit dir.
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON

: "${BASE_START:?set BASE_START=YYYY-MM-DD}"
: "${BASE_END:?set BASE_END=YYYY-MM-DD}"
NJOBS="${SLURM_CPUS_PER_TASK:-32}"

"${PYTHON}" "${HERE}/02_t_climatology_6h.py" "${BASE_START}" "${BASE_END}" "${NJOBS}"
