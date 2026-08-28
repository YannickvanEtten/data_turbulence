#!/bin/bash
# ===========================================================================
# ada/progress.sh -- one line per diagnostics job log.
#
#     ada/progress.sh                 # snapshot
#     watch -n 30 ada/progress.sh     # live
#     ada/progress.sh diag-glob       # only logs whose name matches
#     ada/progress.sh '' --all        # include stale logs from past runs
#
# `tail -f` is the right tool for ONE job. Once an array is running there are
# a dozen logs and the question changes from "what is this job printing" to
# "how far along is each of them". This answers the second one.
#
# Progress comes from the LOGS -- squeue says RUNNING, which after twenty
# minutes is not the information you want. But liveness comes from SQUEUE: a
# log file alone cannot tell you whether its job is still alive, so a cancelled
# or long-finished run would otherwise sit at "loading" forever. Logs whose
# job id is no longer queued and which never reached PEAK RSS are marked
# 'stale' and hidden unless you pass --all.
# ===========================================================================
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
FILTER="${1:-diag-}"
SHOW_ALL="${2:-}"

shopt -s nullglob
logs=(logs/*"$FILTER"*.out)
if [ ${#logs[@]} -eq 0 ]; then
    echo "no logs matching 'logs/*${FILTER}*.out'"
    exit 0
fi

# Job ids SLURM currently knows about (running or pending).
LIVE=" "
if command -v squeue >/dev/null; then
    LIVE=" $(squeue --me -h -o '%i' 2>/dev/null | tr '\n' ' ')"
fi

printf "%-26s %-10s %-10s %s\n" "JOB" "PROGRESS" "RSS" "STATE"
printf "%-26s %-10s %-10s %s\n" "-------------------------" "--------" "--------" "-----"

done_n=0; run_n=0; fail_n=0; stale_n=0

for f in "${logs[@]}"; do
    name=$(basename "$f" .out)
    jobid=${name##*-}                       # diag-glob-1092536_1 -> 1092536_1

    chunk=$(grep -oE 'chunk [0-9]+/[0-9]+' "$f" | tail -1)
    rss=$(grep -oE 'rss [0-9.]+ GB|PEAK RSS  : [0-9.]+ GB' "$f" | tail -1 \
          | grep -oE '[0-9.]+ GB')
    finished=$(grep -c "PEAK RSS" "$f")

    if [ "$finished" -gt 0 ]; then
        vars=$(grep -oE 'variables : [0-9]+ of 21' "$f" | tail -1)
        wall=$(grep -oE 'wall time : .*' "$f" | tail -1 | sed 's/wall time : //')
        if grep -qE '^!!|FAILED:|NO finite values' "$f"; then
            state="DONE - CHECK: $(grep -cE 'FAILED:|NO finite values' "$f") issue(s)"
            fail_n=$((fail_n + 1))
        else
            state="done  ${vars#variables : }  ${wall}"
            done_n=$((done_n + 1))
        fi
        chunk="${chunk:-complete}"

    elif [[ "$LIVE" == *" $jobid "* ]]; then
        # Alive. No chunk line yet means Python is still decoding the GRIB.
        if [ -z "$chunk" ]; then state="loading"; chunk="-"; else state="running"; fi
        run_n=$((run_n + 1))

    else
        # Not finished, not in the queue: cancelled, killed, or a log from an
        # older job whose format this script does not read.
        state="stale (job not in queue)"
        chunk="${chunk:--}"
        stale_n=$((stale_n + 1))
        [ "$SHOW_ALL" != "--all" ] && continue
    fi

    printf "%-26s %-10s %-10s %s\n" "$name" "${chunk#chunk }" "${rss:--}" "$state"
done

echo
line="$done_n done, $run_n in progress, $fail_n need a look"
[ "$stale_n" -gt 0 ] && [ "$SHOW_ALL" != "--all" ] \
    && line="$line, $stale_n stale hidden (--all to show)"
echo "$line"

if command -v squeue >/dev/null; then
    pending=$(squeue --me -h -t PD -o "%i" 2>/dev/null | tr '\n' ' ')
    [ -n "$pending" ] && echo "queued (no log yet): $pending"
fi
