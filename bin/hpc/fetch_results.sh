#!/usr/bin/env bash
# RUN ON THE MAC. Copies the small outputs back (summaries, figures, tables, logs, report),
# NOT the big .h5ad matrices.
#   bash bin/hpc/fetch_results.sh                 # into ~/cln3_msi_outputs_from_ada
#   bash bin/hpc/fetch_results.sh ~/somewhere
set -euo pipefail
source "$(dirname "$0")/_common.sh"
DEST="${1:-$HOME/cln3_msi_outputs_from_ada}"
mkdir -p "$DEST"
rsync -az -e "ssh ${SSH_OPTS[*]}" \
  --include "*/" --exclude "processed/**" --exclude "cache/**" \
  "$REMOTE:$HPC_WORK/outputs/" "$DEST/"
rsync -az -e "ssh ${SSH_OPTS[*]}" "$REMOTE:$HPC_WORK/slurm_logs/" "$DEST/slurm_logs/"
echo "Copied to $DEST"
echo "  summaries: $DEST/summaries   report: $DEST/report/report.html"
