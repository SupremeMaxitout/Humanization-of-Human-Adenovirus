#!/bin/bash
#
# concat_trajectories.sh — join -noappend part files into one trajectory.
#
# Production runs that span several 24 h windows write prod.xtc,
# prod.part0002.xtc, prod.part0003.xtc ... Analysis MUST use the joined
# trajectory, or you will silently analyse only the first segment of one system
# against the full 100 ns of another — which quietly corrupts exactly the RMSF
# comparison the project rests on.
#
# Produces prod_full.xtc per run. Safe to re-run.
#
set -uo pipefail
BASE="${BASE_DIR:-$HOME/md_runs}"

for d in "$BASE"/*/; do
    [ -d "$d" ] || continue
    S=$(basename "$d")
    cd "$d" || continue

    parts=$(ls prod.xtc prod.part*.xtc 2>/dev/null | wc -l)
    if [ "$parts" -eq 0 ]; then
        echo "[skip] $S — no trajectory"; continue
    fi
    if [ -s prod_full.xtc ] && [ prod_full.xtc -nt "$(ls -t prod.part*.xtc prod.xtc 2>/dev/null | head -1)" ]; then
        echo "[done] $S — prod_full.xtc already current"; continue
    fi

    echo "[$S] joining $parts segment(s)"
    # 'c' = continue: take each file's start time from the end of the previous
    yes c | gmx trjcat -f prod.xtc prod.part*.xtc -o prod_full.xtc -settime >/dev/null 2>&1

    if [ -s prod_full.xtc ]; then
        gmx check -f prod_full.xtc 2>&1 | grep -E "^Step" | sed 's/^/    /'
    else
        echo "    [ERROR] trjcat produced nothing"
    fi
done

echo ""
echo "Analysis uses prod_full.xtc. Note that a walltime kill can leave a small"
echo "gap between segments (the run continues past the last written frame), so"
echo "frame spacing may not be perfectly uniform. This does not affect RMSD,"
echo "RMSF or Rg, which are computed per frame."
