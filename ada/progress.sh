#!/bin/bash
# ===========================================================================
# ada/progress.sh -- one line per running/finished diagnostics job.
#
#     ada/progress.sh                 # snapshot
#     watch -n 30 ada/progress.sh     # live
#     ada/progress.sh diag-glob       # only jobs whose log name matches
#
# `tail -f` is the right tool for ONE job. Once an array is running there are
# a dozen logs and the question changes from "what is this job printing" to
# "how far along is each of them". This answers the second one.
#
# Reads the logs, not squeue: the per-chunk progress lines are the only place
# that says how far into a file a task actually is. squeue tells you a task is
# RUNNING, which after twenty minutes is not the information you want.
# ===========================================================================
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
FILTER="${1:-diag-}"

shopt -s nullglob
logs=(logs/*"$FILTER"*.out)
if [ ${#logs[@]} -eq 0 ]; then
    echo "no logs matching 'logs/*${FILTER}*.out'"
    exit 0
fi

printf "%-26s %-12s %-11s %s\n" "JOB" "PROGRESS" "RSS" "STATE"
printf "%-26s %-12s %-11s %s\n" "--------------------------" "------------" "-----------" "-----"

done_n=0; run_n=0; fail_n=0

for f in "${logs[@]}"; do
    name=$(basename "$f" .out)

    # Last "chunk i/n" line, if any. Files run with --chunk-days 0 have a
    # single chunk and will sit at 1/1 for their whole run -- expected.
    chunk=$(grep -oE 'chunk [0-9]+/[0-9]+' "$f" | tail -1)
    rss=$(grep -oE 'rss [0-9.]+ GB|PEAK RSS  : [0-9.]+ GB' "$f" | tail -1 \
          | grep -oE '[0-9.]+ GB')

    if grep -q "PEAK RSS" "$f"; then
        vars=$(grep -oE 'variables : [0-9]+ of 21' "$f" | tail -1 | grep -oE '[0-9]+ of 21')
        wall=$(grep -oE 'wall time : .*' "$f" | tail -1 | sed 's/wall time : //')
        if grep -qE '^!!|FAILED:|NO finite values' "$f"; then
            state="DONE - CHECK: $(grep -cE 'FAILED:|NO finite values' "$f") issue(s)"
            fail_n=$((fail_n + 1))
        else
            state="done  ${vars}  ${wall}"
            done_n=$((done_n + 1))
        fi
        chunk="${chunk:-complete}"
    elif [ -z "$chunk" ]; then
        # Header printed by the shell, but Python has not reached a chunk yet:
        # still decoding the GRIB.
        state="loading"
        chunk="-"
        run_n=$((run_n + 1))
    else
        state="running"
        run_n=$((run_n + 1))
    fi

    printf "%-26s %-12s %-11s %s\n" "$name" "${chunk#chunk }" "${rss:--}" "$state"
done

echo
echo "$done_n done, $run_n in progress, $fail_n need a look"

# Queued array tasks have no log file yet, so they are invisible above.
if command -v squeue >/dev/null; then
    pending=$(squeue --me -h -t PD -o "%i" 2>/dev/null | tr '\n' ' ')
    [ -n "$pending" ] && echo "queued (no log yet): $pending"
fi
