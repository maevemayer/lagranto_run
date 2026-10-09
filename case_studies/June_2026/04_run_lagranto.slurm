#!/bin/bash
#SBATCH --job-name=lagr6h
#SBATCH --partition=zen16           # [CHECK] spiritx partition (sinfo -s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
# ---------------------------------------------------------------------------
# Stage 05: LAGRANTO 15-day backward trajectories for every 6-hourly start.
# PREREQUISITE (one-time, manual): set  set interval = "6"  in
#   $LAGRANTO/goodies/datelist.sh
# Submit:  sbatch --export=ALL,START=2026-05-10,END=2026-06-01 04_run_lagranto.slurm
# ---------------------------------------------------------------------------
set -euo pipefail
# Slurm runs a spool copy of this script, so BASH_SOURCE is useless here:
# use PIPELINE_DIR (exported by submit_all.sh), else the sbatch submit dir.
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON

: "${START:?set START=YYYY-MM-DD}"
: "${END:?set END=YYYY-MM-DD}"
NJOBS="${SLURM_CPUS_PER_TASK:-32}"

"${PYTHON}" "${HERE}/04_run_lagranto_driver.py" "${START}" "${END}" "${NJOBS}"
