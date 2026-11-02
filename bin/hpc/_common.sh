# Shared bits for the bin/hpc scripts (sourced, not run).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_LOCAL="$(cd "$HERE/../.." && pwd)"
# shellcheck disable=SC1091
source "$HERE/settings.sh"

if [[ "$HPC_USER" == your_ada_username* || "$HPC_HOST" == *your-institution* || "$HPC_WORK" == *your_ada_username* ]]; then
  echo "Please edit bin/hpc/settings.sh first (HPC_USER, HPC_HOST, HPC_WORK)." >&2
  exit 2
fi
if [[ "$HPC_WORK" == *" "* ]]; then
  echo "HPC_WORK must not contain spaces: $HPC_WORK" >&2
  exit 2
fi

REMOTE="$HPC_USER@$HPC_HOST"
REMOTE_REPO="$HPC_WORK/repo"
# One ssh connection reused for everything, so you type your password / 2FA only once
SSH_OPTS=(-o ControlMaster=auto -o "ControlPath=$HOME/.ssh/cm-%r@%h:%p" -o ControlPersist=30m)
mkdir -p "$HOME/.ssh" 2>/dev/null || true
