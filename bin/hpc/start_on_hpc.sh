#!/usr/bin/env bash
# RUN ON THE MAC. Logs in to Ada, sets up the software (first time only), downloads the
# reference data, and submits the whole pipeline as SLURM jobs.
#   bash bin/hpc/start_on_hpc.sh
#   FROM=05 bash bin/hpc/start_on_hpc.sh     # resubmit from phase 05 on (earlier outputs reused)
# The jobs keep running on Ada after you close the terminal.
set -euo pipefail
source "$(dirname "$0")/_common.sh"
ssh -t "${SSH_OPTS[@]}" "$REMOTE" "cd '$REMOTE_REPO' && bash bin/hpc/hpc_setup.sh && FROM='${FROM:-01}' bash bin/hpc/submit_all.sh"
