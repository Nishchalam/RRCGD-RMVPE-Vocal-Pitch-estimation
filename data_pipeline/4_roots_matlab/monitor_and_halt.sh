#!/usr/bin/env bash
#
# monitor_and_halt.sh
#
# Watches 1-min system load average. If it stays above LOAD_THRESHOLD for
# SUSTAIN_CHECKS consecutive polls, kills ONLY the MATLAB root-extraction
# processes (the MATLAB binary itself, matched by executable path + its
# distcomp parpool workers) to protect the box from the other running
# GPU-job screens. It deliberately does NOT match on the string
# "run_experiment" alone, since GPU-job scripts also use that filename.
#
# Usage:
#   ./monitor_and_halt.sh [LOAD_THRESHOLD] [POLL_SECONDS] [SUSTAIN_CHECKS]
#
# Defaults: LOAD_THRESHOLD=20, POLL_SECONDS=30, SUSTAIN_CHECKS=3
#   (so it only kills after ~90s of sustained overload, not one noisy spike)

LOAD_THRESHOLD="${1:-20}"
POLL_SECONDS="${2:-30}"
SUSTAIN_CHECKS="${3:-3}"

MATLAB_BIN="/home/Apps/matlab/bin/glnxa64/MATLAB"
LOG_FILE="$(dirname "$0")/monitor_and_halt.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

matlab_pids() {
    # Require BOTH: (1) the executable is the actual MATLAB binary, not
    # any python/GPU process, MATLABWindow helper, or the VSCode MATLAB
    # language-server process that also runs under this same binary, and
    # (2) the command line shows it's the root-extraction job specifically
    # (the "run_experiment" batch driver or its distcomp parpool workers).
    pgrep -f "^${MATLAB_BIN} .*(run_experiment|distcomp_evaluate_filetask)"
}

kill_job() {
    local pids
    pids=$(matlab_pids)
    log "Halting MATLAB process(es): $pids"
    kill -TERM $pids 2>/dev/null
    sleep 5
    pids=$(matlab_pids)
    if [ -n "$pids" ]; then
        log "Process(es) still alive after SIGTERM, sending SIGKILL: $pids"
        kill -KILL $pids 2>/dev/null
    fi
    log "Job halted."
}

log "Starting monitor: threshold=${LOAD_THRESHOLD} poll=${POLL_SECONDS}s sustain=${SUSTAIN_CHECKS} checks (target: $MATLAB_BIN only)"

overload_count=0

while true; do
    if [ -z "$(matlab_pids)" ]; then
        log "Target MATLAB job no longer running. Monitor exiting."
        exit 0
    fi

    load1=$(awk '{print $1}' /proc/loadavg)
    over=$(awk -v l="$load1" -v t="$LOAD_THRESHOLD" 'BEGIN{print (l > t) ? 1 : 0}')

    if [ "$over" -eq 1 ]; then
        overload_count=$((overload_count + 1))
        log "Load $load1 > threshold $LOAD_THRESHOLD (streak: $overload_count/$SUSTAIN_CHECKS)"
    else
        if [ "$overload_count" -gt 0 ]; then
            log "Load back to normal ($load1). Resetting streak."
        fi
        overload_count=0
    fi

    if [ "$overload_count" -ge "$SUSTAIN_CHECKS" ]; then
        log "Sustained overload detected (load=$load1). Triggering halt."
        kill_job
        exit 1
    fi

    sleep "$POLL_SECONDS"
done
