#!/bin/bash
#
# generate_all_jobs.sh — create PBS scripts for every system in systems.conf
#
# Run ONCE from the login node. Safe to re-run: it rewrites job scripts but
# never deletes simulation data.
#
# Environment:
#   FORCE_FIELD  charmm27 (default, what the reference runs used) or charmm36m
#   REPLICATES   number of independent runs per system (default 1)
#   TRIMER_DIR   where the input PDBs live
#   BASE_DIR     where runs are created
#
# On replicates: a single 100 ns run per system is a SCREEN, not converged
# sampling. Differences of the order of run-to-run noise are not real
# differences. Set REPLICATES=3 to generate independent runs (different
# velocity seeds) so you can put error bars on the comparison. That is what a
# reviewer will ask for. It also triples the compute.
#
set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:-$HOME/md_pipeline}"
BASE_DIR="${BASE_DIR:-$HOME/md_runs}"
TRIMER_DIR="${TRIMER_DIR:-$HOME/trimers}"
CONF="${CONF:-$PIPELINE_DIR/systems.conf}"
FORCE_FIELD="${FORCE_FIELD:-charmm27}"
REPLICATES="${REPLICATES:-1}"

[ -s "$CONF" ] || { echo "[FATAL] missing $CONF"; exit 1; }
mkdir -p "$BASE_DIR"

echo "================================================================"
echo "  Generating MD jobs"
echo "  config      : $CONF"
echo "  force field : $FORCE_FIELD"
echo "  replicates  : $REPLICATES"
echo "================================================================"

nok=0; nbad=0
while read -r SYSNAME PDBNAME; do
    [[ -z "${SYSNAME:-}" || "$SYSNAME" == \#* ]] && continue

    INPUT_PDB="$TRIMER_DIR/$PDBNAME"
    if [ ! -s "$INPUT_PDB" ]; then
        echo "  [MISSING] $SYSNAME -> $INPUT_PDB"
        nbad=$((nbad+1)); continue
    fi

    for r in $(seq 1 "$REPLICATES"); do
        if [ "$REPLICATES" -gt 1 ]; then
            RUNNAME="${SYSNAME}_rep${r}"
        else
            RUNNAME="$SYSNAME"
        fi

        WORK="$BASE_DIR/$RUNNAME"
        mkdir -p "$WORK"
        cp -f "$PIPELINE_DIR/mdp/"*.mdp "$WORK/"

        # replicates need independent velocities
        if [ "$REPLICATES" -gt 1 ]; then
            sed -i -e "s/^gen-vel.*/gen-vel = yes/" \
                   -e "s/^gen-seed.*/gen-seed = $((RANDOM + r))/" \
                   "$WORK/prod.mdp"
            grep -q "^gen-temp" "$WORK/prod.mdp" || echo "gen-temp = 310" >> "$WORK/prod.mdp"
        fi

        for t in setup prod; do
            sed -e "s|SYSNAME|$RUNNAME|g" \
                -e "s|INPUT_PDB|$INPUT_PDB|g" \
                -e "s|BASEDIR|$BASE_DIR|g" \
                -e "s|PIPELINEDIR|$PIPELINE_DIR|g" \
                -e "s|FFNAME|$FORCE_FIELD|g" \
                "$PIPELINE_DIR/job_${t}.pbs.template" > "$WORK/job_${t}.pbs"
        done

        chains=$(grep "^ATOM" "$INPUT_PDB" | cut -c22 | sort -u | tr -d '\n')
        echo "  [OK] $RUNNAME  (chains: $chains)"
        nok=$((nok+1))
    done
done < "$CONF"

echo ""
echo "  $nok run(s) generated, $nbad skipped."
echo ""
echo "  NEXT:"
echo "    bash $PIPELINE_DIR/tools/check_status.sh"
echo "    for d in $BASE_DIR/*/; do qsub \"\$d/job_setup.pbs\"; done"
echo "================================================================"
