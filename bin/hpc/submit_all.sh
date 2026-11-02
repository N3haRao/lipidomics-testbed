#!/usr/bin/env bash
# RUNS ON ADA (start_on_hpc.sh calls it). Submits the pipeline as a chain of SLURM jobs:
#
#   01 inventory ─┬─ 02 QC slide E
#                 ├─ lip_neg: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┐
#                 ├─ lip_pos: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┼─ 10 report
#                 └─ met_neg: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┘
#
# Each job starts only if the one before it finished OK (afterok), so a failure stops its
# own chain and the report, and leaves the other modes running.
#   FROM=05 bash bin/hpc/submit_all.sh    # skip 01-04 (their outputs are reused)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/settings.sh"
HPC_WORK="${HPC_WORK:-$HOME/cln3_msi}"   # empty in settings.sh = ~/cln3_msi
export CLN3_ENV_FILE="$HPC_WORK/cln3_env.sh"
[[ -f "$CLN3_ENV_FILE" ]] || { echo "run bin/hpc/hpc_setup.sh first"; exit 2; }
command -v sbatch >/dev/null || { echo "sbatch not found: is Ada a SLURM cluster?"; exit 2; }
FROM="${FROM:-01}"
LOGS="$HPC_WORK/slurm_logs"
mkdir -p "$LOGS"

COMMON=(--parsable --export=ALL)
[[ -n "$SLURM_PARTITION" ]] && COMMON+=(--partition="$SLURM_PARTITION")
[[ -n "$SLURM_ACCOUNT" ]] && COMMON+=(--account="$SLURM_ACCOUNT")
[[ -n "$SLURM_QOS" ]] && COMMON+=(--qos="$SLURM_QOS")

# time limit per phase (edit if a job gets killed for running out of time)
declare -A TIME=([01]=01:00:00 [02]=01:00:00 [03]=08:00:00 [04]=06:00:00 [05]=00:30:00
                 [06]=01:00:00 [07]=02:00:00 [08]=02:00:00 [09]=02:00:00 [10]=00:30:00)
ORDER=(03 04 05 06 07 08 09)

# phase string comparison: run phase P if P >= FROM
want() { [[ "$1" > "$FROM" || "$1" == "$FROM" ]]; }

submit() {   # submit <phase> <mode or ""> <dependency job ids, colon separated or "">
  local phase=$1 mode=$2 dep=$3 extra=()
  [[ -n "$dep" ]] && extra+=(--dependency="afterok:$dep" --kill-on-invalid-dep=yes)
  PHASE="$phase" MODE="$mode" sbatch "${COMMON[@]}" ${extra[@]+"${extra[@]}"} \
    --job-name="cln3_${phase}${mode:+_$mode}" --time="${TIME[$phase]}" \
    --output="$LOGS/%x_%j.out" "$HERE/job_phase.sbatch"
}

dep0=""
if want 01; then
  j=$(submit 01 "" ""); echo "01 inventory          job $j"; dep0=$j
fi
if want 02; then
  j=$(submit 02 "" "$dep0"); echo "02 QC slide E         job $j"
fi
last_ids=()
for mode in $MODES; do
  dep="$dep0"
  for p in "${ORDER[@]}"; do
    want "$p" || continue
    j=$(submit "$p" "$mode" "$dep"); echo "$p $mode              job $j"
    dep=$j
  done
  [[ -n "$dep" ]] && last_ids+=("$dep")
done
deps=$(IFS=:; echo "${last_ids[*]:-}")
j=$(submit 10 "" "$deps"); echo "10 report             job $j"

echo
echo "Submitted. Watch with:   squeue -u $HPC_USER"
echo "Job logs:                $LOGS"
echo "Phase summaries:         $HPC_WORK/outputs/summaries"
echo "From the Mac:            bash bin/hpc/status.sh   /   bash bin/hpc/fetch_results.sh"
