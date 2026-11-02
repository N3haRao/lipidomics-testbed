#!/usr/bin/env bash
# RUN ON THE MAC. Copies the code and the imzML data to Ada.
#   bash bin/hpc/upload_to_hpc.sh            # code + data
#   bash bin/hpc/upload_to_hpc.sh --code     # code only (after I change something)
#
# Uses rsync, so it can be stopped and re-run: it picks up where it left off and skips files
# that are already there and identical. ~40 GB, so the first run can take a while; the
# `caffeinate` keeps the Mac awake meanwhile.
#
# Only "imzML files for biostats" is uploaded. The SCiLS projects, Pilar's results and the
# 09SEP26 cerebellum folder are not needed for the analysis.
set -euo pipefail
source "$(dirname "$0")/_common.sh"

CODE_ONLY=0
[[ "${1:-}" == "--code" ]] && CODE_ONLY=1
SRC_DATA="$LOCAL_DATA_ROOT/imzML files for biostats"

if [[ $CODE_ONLY -eq 0 ]]; then
  if [[ ! -d "$SRC_DATA" ]]; then
    echo "Not found: $SRC_DATA"
    echo "Fix LOCAL_DATA_ROOT in bin/hpc/settings.sh"
    exit 2
  fi
  n=$(find "$SRC_DATA" -type f \( -name "*.imzML" -o -name "*.ibd" \) | wc -l | tr -d ' ')
  echo "Local data: $n imzML/ibd files (expected 36), $(du -sh "$SRC_DATA" | cut -f1)"
  [[ "$n" -eq 36 ]] || echo "WARNING: expected 36 files, found $n. Continuing anyway."
fi

echo "Connecting to $REMOTE (you may be asked for your password / 2FA once)"
ssh "${SSH_OPTS[@]}" "$REMOTE" "mkdir -p '$REMOTE_REPO' '$HPC_WORK/data' '$HPC_WORK/slurm_logs'"

echo "== uploading code to $REMOTE:$REMOTE_REPO"
rsync -az --delete -e "ssh ${SSH_OPTS[*]}" \
  --exclude ".git/" --exclude "synthetic_data/" --exclude "outputs/" --exclude "__pycache__/" \
  --exclude "*.egg-info/" --exclude ".ipynb_checkpoints/" --exclude ".DS_Store" \
  "$REPO_LOCAL/" "$REMOTE:$REMOTE_REPO/"

if [[ $CODE_ONLY -eq 0 ]]; then
  echo "== uploading data to $REMOTE:$HPC_WORK/data/ (this is the long part)"
  # The destination has no spaces on purpose; the folder names with spaces are copied as they are.
  caffeinate -i rsync -a --partial --progress -e "ssh ${SSH_OPTS[*]}" \
    "$SRC_DATA" "$REMOTE:$HPC_WORK/data/"
  echo "== checking what arrived"
  ssh "${SSH_OPTS[@]}" "$REMOTE" "cd '$HPC_WORK/data' && find . -type f \( -name '*.imzML' -o -name '*.ibd' \) | wc -l && du -sh ."
  echo "(first number should be 36)"
fi
echo
echo "Upload done. Next, from the Mac:  bash bin/hpc/start_on_hpc.sh"
