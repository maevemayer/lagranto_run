#!/bin/bash
#SBATCH --job-name=sfiles6h
#SBATCH --partition=zen16           # [CHECK] spiritx partition (sinfo -s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/homedata/%u/lagranto_run/slurm_logs/%x-%j.out
# ---------------------------------------------------------------------------
# Stage 03: build S-files, then append hyam/hybm (stage 04).
# Submit:  sbatch --export=ALL,START=2026-05-15,END=2026-06-01 03_build_s_files.slurm
# ---------------------------------------------------------------------------
set -euo pipefail
# Slurm runs a spool copy of this script, so BASH_SOURCE is useless here:
# use PIPELINE_DIR (exported by submit_all.sh), else the sbatch submit dir.
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON

: "${START:?set START=YYYY-MM-DD}"
: "${END:?set END=YYYY-MM-DD}"
NJOBS="${SLURM_CPUS_PER_TASK:-32}"

"${PYTHON}" "${HERE}/03_build_s_files_6h.py" "${START}" "${END}" "${NJOBS}"

# Stage 04 is folded into stage 03: hyam/hybm are written straight into each S
# file by 03_build_s_files_6h.py, and the old +1 h rename is obsolete (6-hourly
# snapshots are true to their 00/06/12/18 timestamps). No NCO pass needed.
