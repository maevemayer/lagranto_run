#!/bin/bash
#SBATCH --job-name=cmip6chk
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/homedata/%u/cmip6_check/%x-%j.out
# ---------------------------------------------------------------------------
# Read-only integrity check of /bdd/CMIP6 6hrLev hus/ta/ua/va.
# Submit:  mkdir -p /homedata/$USER/cmip6_check && sbatch tools/check_cmip6_6hrlev.slurm
# Extra args go to the script, e.g. --exp ssp585 --activity ScenarioMIP
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
"${PYTHON}" "${HERE}/tools/check_cmip6_6hrlev.py" --njobs "${SLURM_CPUS_PER_TASK:-32}" "$@"
