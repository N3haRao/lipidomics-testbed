#!/usr/bin/env bash
# End-to-end test on a small fake dataset (no real data needed, no internet needed).
#   bash bin/run_synthetic_test.sh
# Makes synthetic_data/ (if missing), wipes synthetic_data/outputs, runs every phase with
# config.test.yaml, then checks that the planted effects were found. Takes ~5 min.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python}"
export PYTHON="$PY"

if [[ ! -f synthetic_data/truth.csv || "${FRESH:-0}" == "1" ]]; then
  rm -rf synthetic_data
  "$PY" scripts/make_synthetic_data.py --out synthetic_data
fi
rm -rf synthetic_data/outputs synthetic_data/cache

# 00 with the fake LIPID MAPS / HMDB files instead of downloading
"$PY" scripts/00_setup.py --config config.test.yaml \
  --lipidmaps-file synthetic_data/lmsd_fake.zip --hmdb-file synthetic_data/hmdb_fake.zip
FROM=01 bash bin/run_all.sh config.test.yaml
"$PY" -m pytest -q tests
"$PY" tests/check_synthetic_results.py --config config.test.yaml
