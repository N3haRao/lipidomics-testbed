#!/usr/bin/env bash
# RUN ON THE MAC. Copies the code and the imzML data to Ada.
#   bash bin/hpc/upload_to_hpc.sh            # code + data
#   bash bin/hpc/upload_to_hpc.sh --code     # code only (after a git pull / config change)
#
# Uses rsync, so it can be stopped (Ctrl+C) and re-run: it picks up where it left off and
# skips files that are already there and identical. The data is ~40 GB; while it copies,
# `caffeinate` keeps the Mac awake (keep it plugged in, lid open).
#
# Only "imzML files for biostats" is uploaded (36 files). The SCiLS projects (.slx/.sbd),
# .srd/.mca files, SCiLS exports, Pilar's results and the 09SEP26 cerebellum folder stay on
# the Mac: the analysis doesn't use them.
#
# Where things go on Ada (HPC_WORK, default ~/cln3_msi = /home/nehar/cln3_msi):
#   ~/cln3_msi/repo/                           the code (mirrors this folder)
#   ~/cln3_msi/data/imzML files for biostats/  the data
set -euo pipefail
source "$(dirname "$0")/_common.sh"

CODE_ONLY=0
[[ "${1:-}" == "--code" ]] && CODE_ONLY=1
SRC_DATA="$LOCAL_DATA_ROOT/imzML files for biostats"
RSYNC_SSH="ssh ${SSH_OPTS[*]}"
KEEP_AWAKE=()
command -v caffeinate >/dev/null 2>&1 && KEEP_AWAKE=(caffeinate -i)

if [[ $CODE_ONLY -eq 0 ]]; then
  if [[ ! -d "$SRC_DATA" ]]; then
    echo "Not found: $SRC_DATA"
    echo "Fix LOCAL_DATA_ROOT in bin/hpc/settings.sh"
    exit 2
  fi
  n=$(find "$SRC_DATA" -type f \( -name "*.imzML" -o -name "*.ibd" \) | wc -l | tr -d ' ')
  local_kb=$(du -sk "$SRC_DATA" | cut -f1)
  echo "Local data: $n imzML/ibd files (expected 36), $(du -sh "$SRC_DATA" | cut -f1)"
  [[ "$n" -eq 36 ]] || echo "WARNING: expected 36 files, found $n. Continuing anyway."
fi

echo "Connecting to $REMOTE"
ssh "${SSH_OPTS[@]}" "$REMOTE" "mkdir -p '$REMOTE_REPO' '$HPC_WORK/data' '$HPC_WORK/slurm_logs' && echo \"free space on Ada: \$(df -h '$HPC_WORK' | awk 'NR==2{print \$4}')\""

echo "== uploading code to $REMOTE:$REMOTE_REPO"
# --delete makes Ada's copy identical to the Mac's (edit code/config on the Mac, not on Ada).
# synthetic_data/ and outputs/ are excluded, so a test run on Ada is never deleted.
rsync -az --delete -e "$RSYNC_SSH" \
  --exclude ".git/" --exclude "synthetic_data/" --exclude "outputs/" --exclude "__pycache__/" \
  --exclude "*.egg-info/" --exclude ".ipynb_checkpoints/" --exclude ".DS_Store" \
  "$REPO_LOCAL/" "$REMOTE:$REMOTE_REPO/"
echo "code uploaded"

if [[ $CODE_ONLY -eq 0 ]]; then
  echo "== uploading data to $REMOTE:$HPC_WORK/data/ (this is the long part)"
  # The destination path has no spaces on purpose; folder names with spaces are copied as is.
  ${KEEP_AWAKE[@]+"${KEEP_AWAKE[@]}"} rsync -a --partial --progress -e "$RSYNC_SSH" \
    "$SRC_DATA" "$REMOTE:$HPC_WORK/data/"

  echo "== checking what arrived on Ada"
  remote=$(ssh "${SSH_OPTS[@]}" "$REMOTE" "cd '$HPC_WORK/data/imzML files for biostats' && \
    echo \$(find . -type f \( -name '*.imzML' -o -name '*.ibd' \) | wc -l) \$(du -sk . | cut -f1)")
  r_n=${remote% *}
  r_kb=${remote#* }
  echo "on Ada: $r_n files, $((r_kb / 1024 / 1024)) GB   |   on the Mac: $n files, $((local_kb / 1024 / 1024)) GB"
  if [[ "$r_n" -eq "$n" ]]; then
    echo "OK: all files arrived (sizes can differ by a few % between Mac and Linux disks, that's normal)"
  else
    echo "Not all files are there yet: run this script again, it continues where it stopped."
    exit 1
  fi
fi
echo
echo "Done. Next (see docs/ADA_RUN_GUIDE.md):  bash bin/hpc/start_on_hpc.sh"
