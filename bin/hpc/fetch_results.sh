#!/usr/bin/env bash
# RUN ON THE MAC. Copies the small outputs back (summaries, figures, tables, logs, report),
# NOT the big .h5ad matrices.
#   bash bin/hpc/fetch_results.sh                 # into ~/cln3_msi_outputs_from_ada
#   bash bin/hpc/fetch_results.sh ~/somewhere
#   RUN=flipJ bash bin/hpc/fetch_results.sh        # the slide J sensitivity run -> ~/cln3_msi_outputs_from_ada_flipJ
set -euo pipefail
source "$(dirname "$0")/_common.sh"
DEST="${1:-$HOME/cln3_msi_outputs_from_ada${RUN:+_$RUN}}"
mkdir -p "$DEST"
rsync -az -e "ssh ${SSH_OPTS[*]}" \
  --include "*/" --exclude "processed/**" --exclude "cache/**" \
  "$REMOTE:$REMOTE_OUT/" "$DEST/"
rsync -az -e "ssh ${SSH_OPTS[*]}" "$REMOTE:$HPC_WORK/slurm_logs/" "$DEST/slurm_logs/"
echo "Copied to $DEST"
echo "  summaries: $DEST/summaries   report: $DEST/report/report.html"
