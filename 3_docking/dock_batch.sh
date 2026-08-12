#!/bin/bash
#
# dock_batch.sh — generate one PBS job per ligand for HADDOCK 2.5 docking
#                 against a single immune-system component.
#
# This is the REFERENCE IMPLEMENTATION: the script that actually ran the
# docking campaign this project is built on, generalised so any target can be
# supplied. Nothing about the immune panel is hard-coded.
#
# Usage — every target is one invocation:
#
#   TARGET_NAME=HD5 \
#   TARGET_PDB=$HOME/immune_system_components/HD5.pdb \
#   TARGET_AIR=$HOME/trimer_air/HD5_trimer_air.tbl \
#   bash dock_batch.sh
#
# Then:
#   cd $BASE_DIR
#   for f in jobs/*.pbs; do qsub "$f"; done
#
# All configuration is by environment variable so the script stays
# dependency-free on a login node (no python, no yq). Defaults in the CONFIG
# block below are overridden by anything exported before the call.
#
# -----------------------------------------------------------------------------
# CRITICAL FIXES carried over from the original run. Each cost a full queue
# cycle to find, and each is preserved deliberately:
#
#   1. The result-collection block does NOT inherit `set -e`. Every cp is
#      independent, and a failed glob must not abort the rest. THIS was the bug
#      that left work_pbs/ full of structures while results/ stayed empty.
#   2. The PBS log directory is created BEFORE submission, not by the job.
#   3. `$HOME` does not expand inside PBS `-o` directives — absolute paths only.
#   4. `set -e` is scoped to setup and docking, never to post-processing.
#   5. The FIRST run_haddock.py call exits non-zero BY DESIGN (it only writes
#      run.cns). Wrapping it in `set +e` is mandatory; treating it as failure
#      aborts every job.
#   6. Jobs are idempotent: a ligand whose water-refined results are already
#      collected is skipped, so resubmitting the whole batch is always safe.
#   7. PBS job names are sanitised and truncated to 15 chars — PBS rejects
#      long or unusual names.
# -----------------------------------------------------------------------------

set -euo pipefail

# ========================= CONFIG (override by env) =========================
TARGET_NAME="${TARGET_NAME:?Set TARGET_NAME, e.g. TARGET_NAME=HD5}"
TARGET_PDB="${TARGET_PDB:?Set TARGET_PDB, e.g. \$HOME/immune_system_components/HD5.pdb}"
TARGET_AIR="${TARGET_AIR:?Set TARGET_AIR, e.g. \$HOME/trimer_air/HD5_trimer_air.tbl}"

# Structures to dock. MUST be single-chain (see note below).
LIGAND_DIR="${LIGAND_DIR:-$HOME/trimers}"

# Where everything for this target lives.
BASE_DIR="${BASE_DIR:-$HOME/dock_${TARGET_NAME}_run}"

HADDOCK_DIR="${HADDOCK_DIR:-$HOME/software/haddock2.5-2025-08}"

# Resources per job
NCPU="${NCPU:-16}"
MEM="${MEM:-64GB}"
WALLTIME="${WALLTIME:-08:00:00}"
QUEUE="${QUEUE:-normal}"

# Sampling. HADDOCK's own recommended production values are 1000/200/200.
# The pilot runs for this project used 50/25/25 — far cheaper, but too few
# models for cluster statistics to mean much. Choose deliberately and say
# which you used.
N_IT0="${N_IT0:-1000}"
N_IT1="${N_IT1:-200}"
N_WATER="${N_WATER:-200}"

# Clustering (FCC)
CLUST_CUTOFF="${CLUST_CUTOFF:-7.5}"
CLUST_SIZE="${CLUST_SIZE:-4}"

# Optional: restrict to models that passed validation (Tool 2 output).
# PASSED_LIST=reports/validation/passed_models.txt
PASSED_LIST="${PASSED_LIST:-}"
# ============================================================================

# ---------------------------- sanity checks ----------------------------
[ -d "$LIGAND_DIR" ]  || { echo "ERROR: ligand dir not found: $LIGAND_DIR"; exit 1; }
[ -f "$TARGET_PDB" ]  || { echo "ERROR: target PDB not found: $TARGET_PDB"; exit 1; }
[ -f "$TARGET_AIR" ]  || { echo "ERROR: AIR file not found: $TARGET_AIR"; exit 1; }
[ -d "$HADDOCK_DIR" ] || { echo "ERROR: HADDOCK not found: $HADDOCK_DIR"; exit 1; }

if ! grep -qi "assign" "$TARGET_AIR"; then
    echo "ERROR: $TARGET_AIR contains no 'assign' restraints."
    echo "       AIR files are an INPUT you supply — they define which residues"
    echo "       should form the interface. See docs/DOCKING.md."
    exit 1
fi

# Fix #2: create log directory up front, not inside the job.
mkdir -p "$BASE_DIR/logs" "$BASE_DIR/results" "$BASE_DIR/work_pbs" "$BASE_DIR/jobs"

# Stage target PDB + AIR locally so each PBS job is self-contained.
cp -f "$TARGET_PDB" "$BASE_DIR/${TARGET_NAME}.pdb"
cp -f "$TARGET_AIR" "$BASE_DIR/${TARGET_NAME}_air.tbl"

cd "$BASE_DIR"

echo "================================================================"
echo "  Target      : $TARGET_NAME"
echo "  Target PDB  : $TARGET_PDB"
echo "  AIR file    : $TARGET_AIR"
echo "  Ligands     : $LIGAND_DIR"
echo "  Base dir    : $BASE_DIR"
echo "  Sampling    : it0=$N_IT0  it1=$N_IT1  water=$N_WATER"
echo "  Clustering  : FCC cutoff=$CLUST_CUTOFF  min size=$CLUST_SIZE"
echo "  Resources   : $NCPU cpus, $MEM, $WALLTIME, queue=$QUEUE"
echo "================================================================"
echo

count=0
skipped=0

for pdb_file in "$LIGAND_DIR"/*.pdb; do
    [ -e "$pdb_file" ] || continue
    pdb_name=$(basename "$pdb_file" .pdb)

    # Optional filter: only dock models that passed Tool 2 validation.
    if [ -n "$PASSED_LIST" ] && [ -f "$PASSED_LIST" ]; then
        if ! grep -qx "$pdb_name" "$PASSED_LIST"; then
            echo "  [skip] $pdb_name — not in $PASSED_LIST"
            skipped=$((skipped + 1))
            continue
        fi
    fi

    # HADDOCK treats each input as ONE rigid body. A three-chain trimer would
    # be docked as three separate molecules, which is not what you want. Use
    # the single-chain form; the multi-chain form is for MD only.
    n_chains=$(grep "^ATOM" "$pdb_file" | cut -c22 | sort -u | tr -d '\n' | wc -c)
    if [ "$n_chains" -gt 1 ]; then
        echo "  [WARN] $pdb_name has $n_chains chains — HADDOCK expects one rigid body."
        echo "         Use the single-chain output from Tool 1 stage 07."
    fi

    # Fix #7: sanitise and truncate the job name.
    job_name=$(echo "${TARGET_NAME}_${pdb_name}" | sed 's/[^a-zA-Z0-9_-]/_/g' | cut -c1-15)

    cat > "jobs/job_${pdb_name}.pbs" << PBSEOF
#!/bin/bash
#PBS -N ${job_name}
#PBS -q ${QUEUE}
#PBS -l select=1:ncpus=${NCPU}:mem=${MEM}
#PBS -l walltime=${WALLTIME}
#PBS -j oe
#PBS -o ${BASE_DIR}/logs/${pdb_name}.log

# =====================================================================
# ENGINE: HADDOCK 2.5  (run.param + run.cns interface)
# Docking: ${TARGET_NAME} <-> ${pdb_name}
#
# Fix #3: absolute path above — \$HOME does not expand in PBS -o.
# Resubmitting this script is safe (see Phase 1).
# =====================================================================

module load python/3.10.9

HADDOCK_DIR="${HADDOCK_DIR}"
BASE_DIR="${BASE_DIR}"
TARGET_PROTEIN="${BASE_DIR}/${TARGET_NAME}.pdb"
AIR_FILE="${BASE_DIR}/${TARGET_NAME}_air.tbl"
PDB_FILE="${pdb_file}"
RESULTS_DIR="${BASE_DIR}/results"
WORK_DIR="${BASE_DIR}/work_pbs"
NCPU="${NCPU}"
PDB_NAME="${pdb_name}"

PYTHON="python3"
RESULT_DIR="\$RESULTS_DIR/\$PDB_NAME"

# ===== Phase 1: idempotency (fix #6) =====
if [ -f "\$RESULT_DIR/it1/file.list" ] && \\
   [ "\$(ls "\$RESULT_DIR/it1/"*.pdb 2>/dev/null | wc -l)" -ge ${N_WATER} ]; then
    echo "[SKIP] \$PDB_NAME already complete"
    exit 0
fi

# ===== Phase 2: setup — strict, fail fast on bad config (fix #4) =====
set -euo pipefail

source "\$HADDOCK_DIR/haddock_configure.sh"
mkdir -p "\$RESULTS_DIR" "\$WORK_DIR"

run_dir="\$WORK_DIR/run_\$PDB_NAME"
mkdir -p "\$run_dir"
cd "\$run_dir"

cat > run.param << PARAM
HADDOCK_DIR=\$HADDOCK_DIR
N_COMP=2
PDB_FILE1=\$PDB_FILE
PDB_FILE2=\$TARGET_PROTEIN
PROJECT_DIR=./docking_\$PDB_NAME
RUN_NUMBER=1
AMBIG_TBL=\$AIR_FILE
STRUCTURES_0=${N_IT0}
STRUCTURES_1=${N_IT1}
WATERREFINE=${N_WATER}
PARAM

# ===== Phase 3: HADDOCK setup pass (fix #5) =====
# This call returns non-zero BY DESIGN; it only generates run.cns.
set +e
\$PYTHON "\$HADDOCK_DIR/haddock/run_haddock.py"
set -e

if [ ! -d "docking_\$PDB_NAME/run1" ]; then
    echo "[ERROR] HADDOCK setup failed for \$PDB_NAME — no run1 directory" >&2
    exit 11
fi

cd "docking_\$PDB_NAME/run1"

# ===== Phase 4: patch run.cns =====
sed -i \\
    -e "s/structures_0=1000/structures_0=${N_IT0}/g" \\
    -e "s/structures_1=200/structures_1=${N_IT1}/g" \\
    -e "s/waterrefine=200/waterrefine=${N_WATER}/g" \\
    -e "s/cpunumber_1=[0-9]*/cpunumber_1=\$NCPU/g" \\
    -e "s/clust_cutoff=.*/clust_cutoff=${CLUST_CUTOFF};/g" \\
    -e "s/clust_size=.*/clust_size=${CLUST_SIZE};/g" \\
    run.cns

# ===== Phase 5: dock =====
echo "[INFO] docking \$PDB_NAME vs ${TARGET_NAME} started \$(date)"
\$PYTHON "\$HADDOCK_DIR/haddock/run_haddock.py"
echo "[INFO] docking finished \$(date)"

# ===== Phase 6: collect — NON-strict ON PURPOSE (fix #1) =====
# Each cp is independent. A missing file.list must NOT prevent the structures
# themselves from being collected. This is the fix for the original bug where
# work_pbs/ had everything and results/ stayed empty.
set +e

mkdir -p "\$RESULT_DIR/it0" "\$RESULT_DIR/it1"

src_it0="structures/it0"
if [ -d "\$src_it0" ]; then
    n_it0=\$(ls "\$src_it0/"*.pdb 2>/dev/null | wc -l)
    [ "\$n_it0" -gt 0 ] && cp "\$src_it0/"*.pdb "\$RESULT_DIR/it0/" 2>/dev/null
    [ -f "\$src_it0/file.list" ] && cp "\$src_it0/file.list" "\$RESULT_DIR/it0/"
    echo "[INFO] collected \$n_it0 it0 structures"
else
    echo "[WARN] it0 directory missing" >&2
fi

src_water="structures/it1/water"
n_water=0
if [ -d "\$src_water" ]; then
    n_water=\$(ls "\$src_water/"*.pdb 2>/dev/null | wc -l)
    [ "\$n_water" -gt 0 ] && cp "\$src_water/"*.pdb "\$RESULT_DIR/it1/" 2>/dev/null
    [ -f "\$src_water/file.list" ] && cp "\$src_water/file.list" "\$RESULT_DIR/it1/"
else
    echo "[ERROR] water refinement directory missing: \$src_water" >&2
fi

# Cluster analysis — required for cluster-level ranking in stage 03.
[ -d "structures/it1/analysis" ] && cp -r "structures/it1/analysis" "\$RESULT_DIR/" 2>/dev/null

# ===== Phase 7: validate =====
if [ "\$n_water" -lt ${N_WATER} ]; then
    echo "[ERROR] expected ${N_WATER} water-refined PDBs, got \$n_water" >&2
    echo "        structures may still be under \$run_dir —" >&2
    echo "        rescue with: python3 3_docking/02_collect_results.py <config>" >&2
    exit 12
fi

if [ ! -f "\$RESULT_DIR/it1/file.list" ]; then
    echo "[ERROR] no file.list in \$RESULT_DIR/it1 — ranking will not work" >&2
    exit 13
fi

echo "[SUCCESS] \$PDB_NAME: \$n_water water-refined structures"
exit 0
PBSEOF

    count=$((count + 1))
done

echo
echo "================================================================"
echo "  Generated $count PBS job(s) in $BASE_DIR/jobs"
[ "$skipped" -gt 0 ] && echo "  Skipped $skipped ligand(s) not in \$PASSED_LIST"
echo "================================================================"
echo
echo "  Submit:"
echo "    for f in $BASE_DIR/jobs/*.pbs; do qsub \"\$f\"; done"
echo
echo "  Monitor:"
echo "    qstat -u \$USER"
echo
echo "  If results/ is empty but work_pbs/ has structures:"
echo "    python3 3_docking/02_collect_results.py <config>"
echo
echo "  Score once all jobs finish:"
echo "    python3 3_docking/03_score_report.py <config>"
