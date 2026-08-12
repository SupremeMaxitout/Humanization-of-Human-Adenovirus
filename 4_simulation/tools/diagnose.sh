#!/bin/bash
#
# diagnose.sh <run> — answer "why did this stop?" in one command.
#
# Ask the questions in this order. Jumping straight to "the structure must be
# bad" is what cost two weeks the first time round.
#
set -uo pipefail
S="${1:?Usage: bash diagnose.sh <run_name>}"
d="${BASE_DIR:-$HOME/md_runs}/$S"
[ -d "$d" ] || { echo "no such run: $d"; exit 1; }

echo "=================================================================="
echo " DIAGNOSIS: $S"
echo "=================================================================="

echo ""
echo "--- 1. STAGE LADDER (first MISSING file is where it stopped) ---"
for f in processed.gro boxed.gro solvated.gro ionised.gro em.gro \
         nvt.gro npt_warmup.gro npt.gro prod.tpr; do
    [ -s "$d/$f" ] && echo "   [have] $f" || echo "   [MISS] $f"
done
ls "$d"/prod*.gro >/dev/null 2>&1 && echo "   [have] prod*.gro (production reached nsteps)" \
                                  || echo "   [MISS] prod*.gro"

echo ""
echo "--- 2. KILLED, OR CRASHED? ---"
if grep -qi "walltime.*exceeded" "$d"/*.log 2>/dev/null; then
    echo "   >> WALLTIME KILL. NOT a physics problem."
    echo "   >> ACTION: resubmit. -cpi resumes from the checkpoint."
elif grep -qi "Run time exceeded" "$d"/*.log 2>/dev/null; then
    echo "   >> -maxh graceful stop. Checkpoint written."
    echo "   >> ACTION: resubmit to continue."
elif ls "$d"/step*.pdb >/dev/null 2>&1; then
    echo "   >> CRASH dumps present:"
    ls -la --time-style=+%Y-%m-%d_%H:%M "$d"/step*.pdb 2>/dev/null | head -3 | sed 's/^/      /'
    echo "   >> CHECK TIMESTAMPS — these go stale and survive from old attempts."
else
    echo "   >> no walltime kill and no crash dumps found"
fi

echo ""
echo "--- 3. ERROR TEXT (grompp errors are usually the real answer) ---"
grep -hiE "Fatal error|Error in user input|cannot|does not exist|not found|Segmentation|Cannot do appending" \
    "$d"/grompp_*.log "$d"/mdrun_*.log "$d"/*_pbs.log 2>/dev/null | sort -u | tail -12

echo ""
echo "--- 4. INSTABILITY SIGNATURE ---"
for L in nvt npt_warmup npt prod; do
    for f in "$d/$L.log" "$d/$L".part*.log; do
        [ -s "$f" ] || continue
        n=$(grep -c "LINCS WARNING" "$f" 2>/dev/null || echo 0)
        s=$(grep -c "can not be settled" "$f" 2>/dev/null || echo 0)
        last=$(grep -E "^ +Step +Time" -A1 "$f" 2>/dev/null | tail -1 | awk '{print $1}')
        echo "   $(basename "$f"): LINCS=$n  settle-fail=$s  last step=${last:-?}"
    done
done

echo ""
echo "--- 5. EXPLODING ATOMS -> RESIDUE ---"
atoms=$(grep -hA6 "bonds that rotated more than" "$d"/npt*.log "$d"/prod*.log 2>/dev/null \
        | grep -oE "^ +[0-9]+ +[0-9]+" | awk '{print $1}' | sort -un | head -6)
if [ -n "$atoms" ] && [ -s "$d/processed.gro" ]; then
    echo "   (atom N is on line N+2 of processed.gro)"
    for a in $atoms; do
        printf "   atom %-7s -> %s\n" "$a" "$(sed -n "$((a+2))p" "$d/processed.gro")"
    done
else
    echo "   none found"
fi

echo ""
echo "=================================================================="
echo " INTERPRETATION  (docs/SIMULATION.md has the full table)"
echo "   walltime / -maxh stop      -> resubmit, nothing is wrong"
echo "   'Cannot do appending'      -> PBS overwrote GROMACS's log."
echo "                                 mv prod.log prod_old.log, add -noappend"
echo "   grompp 'does not exist'    -> a file is missing, not a physics bug"
echo "   LINCS at ~1000-2000 steps of RESTRAINED NPT"
echo "                              -> check refcoord-scaling = com in npt.mdp"
echo "   water 'can not be settled' -> local bad contact; see section 5"
echo "   Epot = inf during EM       -> overlapping atoms in the INPUT structure"
echo "=================================================================="
