#!/bin/bash
#SBATCH --job-name=cmip6chk
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --signal=B:USR1@600
#SBATCH --output=/homedata/%u/cmip6_check/%x-%j.out
# ---------------------------------------------------------------------------
# Read-only integrity check of /bdd/CMIP6 6hrLev hus/ta/ua/va.
# Submit:  mkdir -p /homedata/$USER/cmip6_check && sbatch tools/check_cmip6_6hrlev.slurm
# Default = all activities/experiments. Extra args go to the script, e.g.
#   --quick                              filename inventory only (fast)
#   --best                               longest clean member per model/experiment
#   --activity ScenarioMIP --exp 'ssp*'  one slice per job
#   --fresh                              ignore the checkpoint, check everything again
#   --recheck CESM2,MIROC6               discard and redo the results of these models only
# Default 32 cores (half a zen4 node) starts much sooner than a full node;
# for more:  sbatch --cpus-per-task=64 --mem=160G tools/check_cmip6_6hrlev.slurm ...
#
# Resuming: every finished chunk of 20 files is saved to checkpoint_<tag>.jsonl in the
# output dir. Resubmitting the same command skips what is already done.
# 10 min before the time limit the job resubmits itself (at most
# CMIP6CHK_MAXCHAIN=10 times), so a timeout continues on its own. A crash or
# scancel does not resubmit; sbatch the same command again to resume.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
CHAIN="${CMIP6CHK_CHAIN:-0}"
MAXCHAIN="${CMIP6CHK_MAXCHAIN:-10}"
ARGS=("$@")
# continuations resume: drop --fresh and --recheck (they would discard progress again)
CONT=()
skip=0
for a in ${ARGS[@]+"${ARGS[@]}"}; do
    if (( skip )); then skip=0; continue; fi
    case "$a" in
        --fresh|--recheck=*) ;;
        --recheck) skip=1 ;;
        *) CONT+=("$a") ;;
    esac
done

resubmit() {
    if (( CHAIN < MAXCHAIN )); then
        # keep the resources this job was given (sbatch options override #SBATCH)
        sbatch --dependency=afterany:"${SLURM_JOB_ID}" --chdir="${HERE}" \
               --cpus-per-task="${SLURM_CPUS_PER_TASK:-32}" --mem="${SLURM_MEM_PER_NODE:-98304}M" \
               --export=ALL,CMIP6CHK_CHAIN=$((CHAIN + 1)) \
               "${HERE}/tools/check_cmip6_6hrlev.slurm" ${CONT[@]+"${CONT[@]}"}
        echo "time limit near: continuation $((CHAIN + 1))/${MAXCHAIN} queued, it resumes from the checkpoint"
    else
        echo "time limit near: ${MAXCHAIN} continuations used, sbatch the same command again to resume"
    fi
}
trap resubmit USR1

# python runs in the background so bash can act on USR1 while it works
"${PYTHON}" "${HERE}/tools/check_cmip6_6hrlev.py" --njobs "${SLURM_CPUS_PER_TASK:-32}" ${ARGS[@]+"${ARGS[@]}"} &
PID=$!
wait "${PID}"; rc=$?
while kill -0 "${PID}" 2>/dev/null; do        # wait was interrupted by the signal
    wait "${PID}"; rc=$?
done
exit "${rc}"
