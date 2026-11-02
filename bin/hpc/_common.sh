# Shared bits for the bin/hpc scripts (sourced, not run).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_LOCAL="$(cd "$HERE/../.." && pwd)"
# shellcheck disable=SC1091
source "$HERE/settings.sh"

REMOTE="$HPC_USER@$HPC_HOST"
# One ssh connection reused for everything, so you type your password only once
mkdir -p "$HOME/.ssh" 2>/dev/null || true
SSH_OPTS=(-o ControlMaster=auto -o "ControlPath=$HOME/.ssh/cm-%r@%h:%p" -o ControlPersist=30m)

if [[ -z "$HPC_WORK" ]]; then
  # default: ~/cln3_msi on Ada (ask Ada where home is)
  echo "Connecting to $REMOTE (password prompt, if any, comes now)"
  HPC_WORK="$(ssh "${SSH_OPTS[@]}" "$REMOTE" 'echo $HOME')/cln3_msi"
fi
if [[ "$HPC_WORK" == *" "* ]]; then
  echo "HPC_WORK must not contain spaces: $HPC_WORK" >&2
  exit 2
fi
REMOTE_REPO="$HPC_WORK/repo"
echo "Ada work folder: $HPC_WORK"
