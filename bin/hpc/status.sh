#!/usr/bin/env bash
# RUN ON THE MAC. Shows your jobs on Ada and which phase summaries exist so far.
set -euo pipefail
source "$(dirname "$0")/_common.sh"
ssh "${SSH_OPTS[@]}" "$REMOTE" "squeue -u '$HPC_USER' -o '%.10i %.14j %.8T %.10M %.20R'; echo; echo 'summaries written:'; ls '$HPC_WORK/outputs/summaries' 2>/dev/null || echo '(none yet)'; echo; echo 'failed jobs (if any):'; grep -l -E 'Traceback|Error' '$HPC_WORK'/slurm_logs/*.out 2>/dev/null || echo '(none)'"
