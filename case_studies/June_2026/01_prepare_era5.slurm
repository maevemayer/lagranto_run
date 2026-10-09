#!/bin/bash
#SBATCH --job-name=prep6h
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/homedata/%u/lagranto_run/slurm_logs/%x-%j.out
# ---------------------------------------------------------------------------
# Stage 01: ERA5 -> P-files (6-hourly snapshots).
# Submit:  sbatch --export=ALL,START=2026-05-15,END=2026-06-01 01_prepare_era5.slurm
# ---------------------------------------------------------------------------
set -euo pipefail
# Slurm runs a spool copy of this script, so BASH_SOURCE is useless here:
# use PIPELINE_DIR (exported by submit_all.sh), else the sbatch submit dir.
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON

: "${START:?set START=YYYY-MM-DD}"
: "${END:?set END=YYYY-MM-DD}"
NJOBS="${SLURM_CPUS_PER_TASK:-32}"

"${PYTHON}" "${HERE}/01_prepare_era5_6h.py" "${START}" "${END}" "${NJOBS}"
