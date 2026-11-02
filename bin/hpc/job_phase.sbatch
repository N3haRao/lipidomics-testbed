#!/bin/bash
# One pipeline phase as a SLURM job. submit_all.sh sets PHASE (and MODE) and the time limit.
# The defaults below are what one mode needs (see README): 4 cores, 32 GB.
#SBATCH --job-name=cln3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
set +eu  # a stray warning in the env file (e.g. "module: command not found") must not kill the job
source "${CLN3_ENV_FILE:?CLN3_ENV_FILE not set}"
set -eu
cd "$CLN3_REPO"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
echo "phase ${PHASE} mode ${MODE:-all} on $(hostname), $(date)"
if [[ -n "${MODE:-}" ]]; then
  bash bin/run_phase.sh "$PHASE" config.hpc.yaml --mode "$MODE"
else
  bash bin/run_phase.sh "$PHASE" config.hpc.yaml
fi
echo "finished $(date)"
