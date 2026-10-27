#!/usr/bin/env bash
# Run every phase in order, stop at the first error.
#   bash bin/run_all.sh                    # uses config.yaml
#   bash bin/run_all.sh my_config.yaml
#   FROM=05 bash bin/run_all.sh            # restart from phase 05 (earlier outputs are reused)
# 00 is only the validation here (no downloads); run the --download step once by hand.
set -euo pipefail
cd "$(dirname "$0")/.."
CONFIG="${1:-config.yaml}"
FROM="${FROM:-00}"
PHASES=(00 01 02 03 04 04b 05 06 07 08 09 10)

started=0
t0=$(date +%s)
for p in "${PHASES[@]}"; do
  if [[ "$p" == "$FROM" ]]; then started=1; fi
  if [[ $started -eq 0 ]]; then continue; fi
  bash bin/run_phase.sh "$p" "$CONFIG" || { echo; echo "!!! phase $p failed, stopping. See the newest file in <output_root>/logs/"; exit 1; }
done
echo
echo "all phases done in $(( ($(date +%s) - t0) / 60 )) min. Report: <output_root>/report/report.html"
