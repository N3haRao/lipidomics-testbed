#!/usr/bin/env bash
# RUNS ON ADA (start_on_hpc.sh calls it for you), on the login node, because it needs internet.
#   1. finds conda (module, or one already on your PATH), or installs Miniforge into HPC_WORK
#   2. creates the cln3-msi env in HPC_WORK/envs (first time; reused afterwards)
#   3. runs the unit tests (seconds)
#   4. phase 00: validates config + sample sheet, downloads LIPID MAPS, HMDB, Allen atlas (once)
# Writes HPC_WORK/cln3_env.sh, which every SLURM job sources to get the same environment.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
source "$HERE/settings.sh"
HPC_WORK="${HPC_WORK:-$HOME/cln3_msi}"   # empty in settings.sh = ~/cln3_msi
export CLN3_WORK="$HPC_WORK"
ENV_PREFIX="$HPC_WORK/envs/cln3-msi"
cd "$REPO"

# ---- 1. conda ---------------------------------------------------------------------------
CONDA_SH=""
set +u   # module and conda scripts don't like `set -u`
if [[ -n "$CONDA_MODULE" ]] && command -v module >/dev/null 2>&1; then
  module load "$CONDA_MODULE"
fi
if command -v conda >/dev/null 2>&1; then
  CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
elif [[ -f "$HPC_WORK/miniforge3/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="$HPC_WORK/miniforge3/etc/profile.d/conda.sh"
else
  echo "== no conda found: installing Miniforge into $HPC_WORK/miniforge3"
  curl -fsSL -o "$HPC_WORK/miniforge.sh" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
  bash "$HPC_WORK/miniforge.sh" -b -p "$HPC_WORK/miniforge3"
  rm -f "$HPC_WORK/miniforge.sh"
  CONDA_SH="$HPC_WORK/miniforge3/etc/profile.d/conda.sh"
fi
# shellcheck disable=SC1090
source "$CONDA_SH"
set -u

# ---- 2. environment ---------------------------------------------------------------------
if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
  echo "== creating the env in $ENV_PREFIX (5-15 min, first time only)"
  conda env create -p "$ENV_PREFIX" -f environment.yml
fi
set +u; conda activate "$ENV_PREFIX"; set -u
pip install -q -e .
pip install -q -r requirements-optional.txt || echo "(optional packages failed, that's fine)"

cat > "$HPC_WORK/cln3_env.sh" <<ENVEOF
# written by hpc_setup.sh; sourced by every SLURM job
${CONDA_MODULE:+module load $CONDA_MODULE}
source "$CONDA_SH"
conda activate "$ENV_PREFIX"
export CLN3_WORK="$HPC_WORK"
export CLN3_REPO="$REPO"
export MPLBACKEND=Agg
ENVEOF

# ---- 3. quick tests ---------------------------------------------------------------------
echo "== unit tests"
python -m pytest -q tests

# ---- 4. phase 00 ------------------------------------------------------------------------
if [[ -f "$HPC_WORK/outputs/cache/lipidmaps/lipidmaps.parquet" && -d "$HPC_WORK/outputs/cache/atlas" ]]; then
  echo "== reference data already cached; validating config only"
  python scripts/00_setup.py --config config.hpc.yaml
else
  echo "== phase 00 with downloads (LIPID MAPS, HMDB, Allen atlas)"
  python scripts/00_setup.py --config config.hpc.yaml --download || {
    echo "00_setup failed. If a download was blocked, see the README (manual download + --lipidmaps-file)."
    exit 1; }
fi
grep -A4 -i "warnings" "$HPC_WORK/outputs/summaries/00_setup.md" || true
echo "== setup OK"
