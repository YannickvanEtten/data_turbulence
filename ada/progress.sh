#!/bin/bash
# ===========================================================================
# ada/progress.sh -- one line per job log, for the run you are watching.
#
#     bash ada/progress.sh              # the NEWEST array only  <- the default
#     bash ada/progress.sh diag-trend   # logs whose name matches a substring
#     bash ada/progress.sh --all        # every log on disk, oldest first
#     bash ada/progress.sh '' --stale   # newest array, including dead logs
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
# 'stale' and hidden unless you pass --stale.
#
# WHY THE DEFAULT IS "NEWEST ARRAY" AND NOT "EVERYTHING"
# -----------------------------------------------------
# It used to default to the substring `diag-`, which matches every diagnostics
# log ever written -- six arrays and counting. Bash expands a glob in
# lexicographic order, so the oldest finished run printed FIRST and the run you
# submitted thirty seconds ago printed last, below the fold. The script looked
# stuck on a completed job when it was really showing all of history.
#
# The default now resolves the most recently modified log, takes its array id,
# and shows only that array. That is almost always the thing you just
# submitted. Everything else is still reachable by argument.
# ===========================================================================
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1

ARG1="${1:-}"
ARG2="${2:-}"
SHOW_STALE="no"
case "$ARG2" in --stale|--all) SHOW_STALE="yes" ;; esac

shopt -s nullglob

all_logs=(logs/*.out)
if [ ${#all_logs[@]} -eq 0 ]; then
    echo "no logs in logs/ -- nothing has run yet, or you are in the wrong directory"
    echo "  pwd: $(pwd)"
    exit 0
fi

scope=""
if [ "$ARG1" = "--all" ]; then
    logs=("${all_logs[@]}")
    scope="every log on disk"

elif [ -n "$ARG1" ]; then
    logs=(logs/*"$ARG1"*.out)
    scope="logs matching '$ARG1'"
    if [ ${#logs[@]} -eq 0 ]; then
        echo "no logs matching 'logs/*${ARG1}*.out'"
        exit 0
    fi

else
    # Newest log by mtime -> its array id -> only that array.
    newest=$(ls -t logs/*.out 2>/dev/null | head -1)
    base=$(basename "$newest" .out)          # e.g. diag-trend-1092956_3
    jid=${base##*-}                          # 1092956_3
    arr=${jid%%_*}                           # 1092956
    logs=(logs/*"-${arr}"_*.out logs/*"-${arr}".out)
    [ ${#logs[@]} -eq 0 ] && logs=("$newest")
    scope="job ${arr} (newest); 'bash ada/progress.sh --all' for the rest"
fi

# Sort by version so task 2 comes before task 10, not after it.
mapfile -t logs < <(printf '%s\n' "${logs[@]}" | sort -V)

# Job ids SLURM currently knows about (running or pending).
LIVE=" "
if command -v squeue >/dev/null; then
    LIVE=" $(squeue --me -h -o '%i' 2>/dev/null | tr '\n' ' ')"
fi

echo "showing: $scope"
echo
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
    skipped=$(grep -c "already exists, skipping" "$f")

    if [ "$skipped" -gt 0 ]; then
        # 07 and the download jobs detect their own output and exit in a
        # second. That is a success, not a no-op worth hiding.
        state="skipped (output already present)"
        chunk="-"
        done_n=$((done_n + 1))

    elif [ "$finished" -gt 0 ]; then
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
        [ "$SHOW_STALE" != "yes" ] && continue
    fi

    printf "%-26s %-10s %-10s %s\n" "$name" "${chunk#chunk }" "${rss:--}" "$state"
done

echo
line="$done_n done, $run_n in progress, $fail_n need a look"
[ "$stale_n" -gt 0 ] && [ "$SHOW_STALE" != "yes" ] \
    && line="$line, $stale_n stale hidden (pass '' --stale to show)"
echo "$line"

if command -v squeue >/dev/null; then
    pending=$(squeue --me -h -t PD -o "%i" 2>/dev/null | tr '\n' ' ')
    [ -n "$pending" ] && echo "queued (no log yet): $pending"
    [ -z "$(squeue --me -h -o '%i' 2>/dev/null)" ] && echo "queue empty -- nothing of yours is running"
fi
