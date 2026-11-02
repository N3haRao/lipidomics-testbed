#!/bin/bash
# End-to-end test on fake data, as a SLURM job (~10 min). From the repo folder on Ada:
#   sbatch bin/hpc/job_test.sbatch
# Log: ~/cln3_msi/slurm_logs/test_<jobid>.out ; the last line should say "19/19 checks passed".
#SBATCH --job-name=cln3_test
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=%x_%j.out
set -eo pipefail
WORK="${CLN3_WORK:-$HOME/cln3_msi}"
source "$WORK/cln3_env.sh"
cd "$CLN3_REPO"
echo "test on $(hostname), $(date), python $(command -v python)"
bash bin/run_synthetic_test.sh
