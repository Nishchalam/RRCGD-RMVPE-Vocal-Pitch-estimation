#!/usr/bin/env bash
# Auto-detecting GPU scheduler for a sweep config.
#
# Detects how many GPUs are present and how much memory each has, works
# out how many jobs each GPU can hold concurrently from that, and
# launches one `run_experiments.py --only <experiment>` per screen --
# each such call runs that experiment's FULL alpha sweep sequentially
# and is itself resumable (--skip-existing at the experiment level,
# EXP_RESUME=1 by default inside train_single.py at the iteration
# level -- see INSTALL.md). Experiments are launched in the order the
# config lists them, one at a time, onto whichever GPU has a free slot
# -- never batch-launches off a stale nvidia-smi reading (always
# re-checks immediately before each launch).
#
# Every 10 experiments that complete (test_eval.csv appears on disk --
# counted on-disk, never from a log), re-runs export_validation_csv.py
# and compile_best_results.py so the aggregated CSVs stay fresh through
# a long sweep without you having to remember to do it by hand.
#
# Usage:
#   ./auto_gpu_scheduler.sh <config.yaml> [per_job_mib]
#
# per_job_mib (optional, default 6500): how much GPU memory one job
# needs, used both to size each GPU's capacity (total_mem / per_job_mib)
# and as the free-memory gate before every launch. The project's own
# empirical numbers (see the main README/INSTALL.md): ~3.8-6.5GB/job at
# 256 bins depending on technique, ~4.5-9.5GB/job at 512 bins -- raise
# this for a 512-bin config or any full-range (E1/E2) technique.
set -uo pipefail

CONFIG="${1:?Usage: $0 <config.yaml> [per_job_mib]}"
PER_JOB_MIB="${2:-6500}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$HERE")"
LOGDIR="$HERE/logs"
TICK_S=30
SETTLE_S=60
EXPORT_EVERY=10

mkdir -p "$LOGDIR"
cd "$REPO_ROOT"

log() { echo "[$(date '+%F %T')] $*"; }

# ── GPU detection ─────────────────────────────────────────────────────────
n_gpus() { nvidia-smi --query-gpu=index --format=csv,noheader | wc -l; }

gpu_total_mib() { nvidia-smi -i "$1" --query-gpu=memory.total --format=csv,noheader,nounits; }

gpu_free_mib() {
    read -r used total < <(nvidia-smi -i "$1" --query-gpu=memory.used,memory.total --format=csv,noheader,nounits | tr ',' ' ')
    echo $((total - used))
}

# Jobs this scheduler itself launched that are still running on GPU $1.
running_on_gpu() {
    local g="$1" cnt=0 p
    for p in $(pgrep -f "run_experiments.py --config $CONFIG" 2>/dev/null); do
        if tr '\0' '\n' < "/proc/$p/environ" 2>/dev/null | grep -qx "CUDA_VISIBLE_DEVICES=$g"; then
            cnt=$((cnt + 1))
        fi
    done
    echo "$cnt"
}

N_GPUS=$(n_gpus)
if [ "$N_GPUS" -eq 0 ]; then
    echo "No GPUs detected (nvidia-smi returned nothing). Exiting." >&2
    exit 1
fi

declare -a CAPACITY
for ((g = 0; g < N_GPUS; g++)); do
    total=$(gpu_total_mib "$g")
    cap=$((total / PER_JOB_MIB))
    [ "$cap" -lt 1 ] && cap=1
    CAPACITY[g]=$cap
    log "GPU$g: ${total}MiB total -> capacity ${cap} job(s) @ ${PER_JOB_MIB}MiB each"
done

# ── Job list: every experiment name in the config, in listed order ─────────
mapfile -t EXPERIMENTS < <(python3 -c "
import yaml, sys
c = yaml.safe_load(open(sys.argv[1]))
for e in c['experiments']:
    print(e['name'])
" "$CONFIG")
log "Config $CONFIG: ${#EXPERIMENTS[@]} experiment(s) to run"

already_running() {
    screen -ls 2>/dev/null | grep -q "\.auto_${1}\b"
}

count_test_evals() {
    find runs_seeded -name test_eval.csv 2>/dev/null | wc -l
}

last_export_bucket=$(( $(count_test_evals) / EXPORT_EVERY ))

idx=0
while [ "$idx" -lt "${#EXPERIMENTS[@]}" ]; do
    name="${EXPERIMENTS[$idx]}"

    if already_running "$name"; then
        idx=$((idx + 1))
        continue
    fi

    launched=0
    for ((g = 0; g < N_GPUS; g++)); do
        if [ "$(running_on_gpu "$g")" -lt "${CAPACITY[g]}" ] && [ "$(gpu_free_mib "$g")" -ge "$PER_JOB_MIB" ]; then
            log "LAUNCH $name on GPU$g"
            screen -dmS "auto_${name}" bash -c \
                "CUDA_VISIBLE_DEVICES=$g python3 run_experiments.py --config $CONFIG --only $name --skip-existing >> $LOGDIR/${name}.log 2>&1"
            launched=1
            sleep "$SETTLE_S"
            break
        fi
    done

    if [ "$launched" -eq 1 ]; then
        idx=$((idx + 1))
        continue
    fi

    sleep "$TICK_S"

    # Re-check export bucket on every tick, not just on launch, so a long
    # idle stretch (all GPUs full) still gets its export once jobs finish.
    current=$(count_test_evals)
    bucket=$((current / EXPORT_EVERY))
    if [ "$bucket" -ne "$last_export_bucket" ]; then
        log "=== ${current} test_eval.csv file(s) on disk -- refreshing aggregated CSVs ==="
        python3 export_validation_csv.py
        python3 compile_best_results.py
        last_export_bucket=$bucket
    fi
done

log "All ${#EXPERIMENTS[@]} experiment(s) launched. Waiting for the last ones to finish..."
while true; do
    any_running=0
    for ((g = 0; g < N_GPUS; g++)); do
        [ "$(running_on_gpu "$g")" -gt 0 ] && any_running=1
    done
    [ "$any_running" -eq 0 ] && break
    sleep "$TICK_S"
done

log "=== sweep complete -- final CSV export ==="
python3 export_validation_csv.py
python3 compile_best_results.py
log "Done."
