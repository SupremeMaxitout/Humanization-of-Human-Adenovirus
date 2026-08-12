#!/bin/bash
#
# check_status.sh — where does every run actually stand?
#
# Replaces the naive "[OK]/[FAIL] npt.gro" check, which was actively
# misleading: it reported FAIL for jobs that were merely queued, and gave no
# clue which stage had been reached or how far production had got.
#
# Handles -noappend part files: completion is prod*.gro AND last step == nsteps.
#
set -uo pipefail
BASE="${BASE_DIR:-$HOME/md_runs}"
CONF="${CONF:-$HOME/md_pipeline/systems.conf}"

QUEUED=$(qstat -u "$USER" 2>/dev/null | awk 'NR>5 {print $4}' | tr '\n' ' ')

printf "%-16s %-20s %-10s %s\n" "RUN" "STAGE" "QUEUE" "PROGRESS"
printf "%-16s %-20s %-10s %s\n" "---" "-----" "-----" "--------"

for d in "$BASE"/*/; do
    [ -d "$d" ] || continue
    S=$(basename "$d")

    nsteps=$(grep -E "^[[:space:]]*nsteps" "$d/prod.mdp" 2>/dev/null | head -1 | grep -oE "[0-9]+")
    last=$(grep -hE "^ +Step +Time" -A1 "$d"/prod*.log 2>/dev/null | tail -1 | awk '{print $1}')
    ns=""
    [ -n "${last:-}" ] && ns=$(echo "$last" | awk '{printf "%.1f ns", $1*0.002/1000}')

    if ls "$d"/prod*.gro >/dev/null 2>&1 && [ "${last:-0}" = "${nsteps:-x}" ]; then
        stage="PRODUCTION COMPLETE"
    elif [ -s "$d/prod.cpt" ];                    then stage="production"
    elif ls "$d"/npt.gro "$d"/npt.part*.gro >/dev/null 2>&1; then stage="setup done (ready)"
    elif ls "$d"/npt_warmup*.gro >/dev/null 2>&1; then stage="NPT main"
    elif ls "$d"/nvt*.gro >/dev/null 2>&1;        then stage="NPT warm-up"
    elif [ -s "$d/em.gro" ];                      then stage="NVT"
    elif [ -s "$d/ionised.gro" ];                 then stage="EM"
    elif [ -s "$d/processed.gro" ];               then stage="solvation"
    else                                               stage="not started"
    fi

    q="-"
    echo "$QUEUED" | grep -qE "(setup|prod)_${S}\b" && q="RUNNING/Q"

    printf "%-16s %-20s %-10s %s\n" "$S" "$stage" "$q" "$ns"
done

echo ""
echo "Reading this table:"
echo "  QUEUE=RUNNING/Q -> alive. A missing .gro is NOT a failure yet."
echo "  QUEUE=-         -> nothing queued. If not complete, it stopped:"
echo "                     bash tools/diagnose.sh <run>"
echo "  COMPLETE requires prod*.gro AND last step == nsteps. 'Finished mdrun'"
echo "  alone is not enough — it also prints after a -maxh stop."
