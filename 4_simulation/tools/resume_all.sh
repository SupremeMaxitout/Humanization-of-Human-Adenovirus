#!/bin/bash
#
# resume_all.sh — resubmit every run that is unfinished.
#
# Always safe: setup_md.sh skips finished stages, and all mdrun calls use -cpi,
# so nothing already computed is repeated.
#
#   bash resume_all.sh setup
#   bash resume_all.sh prod
#
set -euo pipefail
MODE="${1:?Usage: bash resume_all.sh <setup|prod>}"
BASE="${BASE_DIR:-$HOME/md_runs}"

cd "$HOME"    # never sit inside a directory you might delete

RUNNING=$(qstat -u "$USER" 2>/dev/null | awk 'NR>5 {print $4}' | tr '\n' ' ')

for d in "$BASE"/*/; do
    [ -d "$d" ] || continue
    S=$(basename "$d")

    if echo "$RUNNING" | grep -qE "(setup|prod)_${S}\b"; then
        echo "[skip] $S — already queued"; continue
    fi

    if [ "$MODE" = "setup" ]; then
        if ls "$d"/npt.gro "$d"/npt.part*.gro >/dev/null 2>&1; then
            echo "[done] $S setup complete"
        else
            echo "[submit] $S setup"; qsub "$d/job_setup.pbs"
        fi
    else
        nsteps=$(grep -E "^[[:space:]]*nsteps" "$d/prod.mdp" 2>/dev/null | head -1 | grep -oE "[0-9]+")
        last=$(grep -hE "^ +Step +Time" -A1 "$d"/prod*.log 2>/dev/null | tail -1 | awk '{print $1}')
        if ls "$d"/prod*.gro >/dev/null 2>&1 && [ "${last:-0}" = "${nsteps:-x}" ]; then
            echo "[done] $S production complete"
        elif ! ls "$d"/npt.gro "$d"/npt.part*.gro >/dev/null 2>&1; then
            echo "[block] $S setup not finished"
        else
            echo "[submit] $S production"; qsub "$d/job_prod.pbs"
        fi
    fi
done
