#!/bin/bash
#SBATCH --job-name=record_hw
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/homedata/%u/lagranto_run/slurm_logs/%x-%j.out
# ---------------------------------------------------------------------------
# Record heatwave days of a CMIP6 model over the paper's target box (find_record_heatwave.py).
# Submit from the repo root:
#   sbatch TXx/find_record_heatwave.slurm                  # MRI-ESM2-0 historical
#   sbatch --export=ALL,MODEL=CNRM-CM6-1 TXx/find_record_heatwave.slurm
# ---------------------------------------------------------------------------
set -u
SRC="${SLURM_SUBMIT_DIR:-$(pwd)}/TXx"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
export CMIP6_POOL="${CMIP6_POOL:-/bdd/CMIP6}"
"${PYTHON}" "${SRC}/find_record_heatwave.py" "${MODEL:-MRI-ESM2-0}" --exp "${EXP:-historical}" \
    --outdir /homedata/${USER}/lagranto_run/TXx
