#!/usr/bin/env bash
# Run one phase:   bash bin/run_phase.sh 04 [config.yaml] [extra args, e.g. --mode lip_neg]
# Phases: 00 01 02 03 04 04b 05 06 07 08 09 10
set -euo pipefail
cd "$(dirname "$0")/.."
PHASE="${1:?usage: run_phase.sh <phase> [config] [args...]}"
CONFIG="${2:-config.yaml}"
shift $(( $# >= 2 ? 2 : 1 ))
PY="${PYTHON:-python}"

case "$PHASE" in
  00)  SCRIPT="scripts/00_setup.py" ;;
  01)  SCRIPT="scripts/01_inventory.py" ;;
  02)  SCRIPT="scripts/02_qc_one_slide.py" ;;
  03)  SCRIPT="scripts/03_qc_all.py" ;;
  04)  SCRIPT="scripts/04_peaks_matrix.py" ;;
  04b) SCRIPT="scripts/04b_scils_crosscheck.py" ;;
  05)  SCRIPT="scripts/05_annotate.py" ;;
  06)  SCRIPT="scripts/06_normalize.py" ;;
  07)  SCRIPT="scripts/07_register.py" ;;
  08)  SCRIPT="scripts/08_explore.py" ;;
  09)  SCRIPT="scripts/09_differential.py" ;;
  10)  SCRIPT="scripts/10_report.py" ;;
  *)   echo "unknown phase $PHASE"; exit 2 ;;
esac
echo ">>> phase $PHASE: $SCRIPT --config $CONFIG $*"
"$PY" "$SCRIPT" --config "$CONFIG" "$@"
