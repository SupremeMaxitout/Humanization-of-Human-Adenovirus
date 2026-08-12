#!/bin/bash
#
# setup_md.sh — prepare ONE system for production MD.
#
#   EM -> NVT (100 ps) -> NPT warm-up (50 ps @ 0.5 fs) -> NPT (1 ns)
#
# Usage:
#   bash setup_md.sh <input.pdb> <system_name>
#
# Configuration is by environment variable (see CONFIG block). Defaults match
# the runs this project actually performed.
#
# -----------------------------------------------------------------------------
# DESIGN NOTES — each is a bug that cost real time. See docs/LESSONS.md.
#
#   * PREFLIGHT before anything else. Most historical failures were a missing
#     file or a login-node execution, not physics.
#   * IDEMPOTENT AND RESUMABLE. Every stage is skipped if its .gro exists, and
#     every mdrun uses -cpi. Re-running after any interruption is always safe.
#   * pdb2gmx uses "-merge no". A homotrimer must stay THREE chains; merging
#     creates artificial peptide bonds across protomer junctions that stiffen
#     the very interface an HVR flexibility study measures.
#   * NO conjugate-gradient stage and NO external "pre-relaxation" of the input.
#     Model coordinates go in as they come — vacuum minimisation of a predicted
#     structure manufactures the defects it is meant to remove.
#   * The refcoord-scaling guard refuses to start if a restrained NPT stage has
#     a barostat without "refcoord-scaling = com".
# -----------------------------------------------------------------------------

set -euo pipefail

INPUT_PDB="${1:?Usage: bash setup_md.sh <input.pdb> <system_name>}"
SYSNAME="${2:?Usage: bash setup_md.sh <input.pdb> <system_name>}"

# ========================= CONFIG (override by env) =========================
BASE_DIR="${BASE_DIR:-$HOME/md_runs}"
MDP_DIR="${MDP_DIR:-$HOME/md_pipeline/mdp}"

# FORCE FIELD — your choice, but be accurate about it in your methods.
#   charmm27  ships with GROMACS. This is what the reference runs used.
#   charmm36m better for long disordered loops (i.e. HVRs), but is a separate
#             port from the MacKerell lab that you must install into $GMXLIB.
# Whichever you pick, use the SAME one for every system in a comparison.
FORCE_FIELD="${FORCE_FIELD:-charmm27}"
WATER_MODEL="${WATER_MODEL:-tip3p}"

BOX_PADDING="${BOX_PADDING:-1.5}"   # nm from protein to box edge
SALT_CONC="${SALT_CONC:-0.15}"      # 150 mM NaCl (extracellular)
NTOMP="${NTOMP:-16}"

# GPU offload for the dynamics stages. Set USE_GPU=0 when DEBUGGING a blow-up:
# CPU LINCS prints exactly which atoms and bonds failed, where a GPU run tends
# to segfault opaquely.
USE_GPU="${USE_GPU:-1}"
if [ "$USE_GPU" = "1" ]; then
    ACCEL="-nb gpu -pme gpu -bonded gpu"
else
    ACCEL="-nb cpu -pme cpu"
fi

# Stop gracefully this many hours in and checkpoint. Must be below the PBS
# walltime so the job can exit cleanly rather than being killed.
MAXH="${MAXH:-23.0}"
# ============================================================================

WORK="$BASE_DIR/$SYSNAME"

# ============================ PREFLIGHT ============================
echo "================================================================"
echo "  PREFLIGHT — $SYSNAME"
echo "================================================================"

if hostname | grep -qi "login"; then
    echo "  [FATAL] You are on a LOGIN NODE ($(hostname))."
    echo "          Submit with qsub, or take an interactive compute node:"
    echo "          qsub -I -l select=1:ncpus=16:ngpus=1 -l walltime=04:00:00 -q normal"
    exit 10
fi
echo "  [ok] compute node $(hostname)"

command -v gmx >/dev/null 2>&1 || {
    echo "  [FATAL] gmx not found. module load gromacs/2023.2-gpu"; exit 11; }
echo "  [ok] gmx present"

[ -s "$INPUT_PDB" ] || { echo "  [FATAL] input missing/empty: $INPUT_PDB"; exit 12; }
natoms=$(grep -c "^ATOM" "$INPUT_PDB" || true)
chains=$(grep "^ATOM" "$INPUT_PDB" | cut -c22 | sort -u | tr -d '\n')
echo "  [ok] input $INPUT_PDB ($natoms atoms, chains: $chains)"
if [ "${#chains}" -lt 3 ]; then
    echo "  [WARN] fewer than 3 chains. A hexon homotrimer should have three."
    echo "         A merged single chain models the protomer interface wrongly."
fi

missing=0
for m in em.mdp nvt.mdp npt_warmup.mdp npt.mdp prod.mdp; do
    [ -s "$MDP_DIR/$m" ] || { echo "  [FATAL] missing mdp: $MDP_DIR/$m"; missing=1; }
done
[ "$missing" = "0" ] || exit 13
echo "  [ok] all 5 mdp files present"

# THE guard. A missing refcoord-scaling with -DPOSRES + barostat is the single
# most expensive bug in this project's history.
for m in npt_warmup.mdp npt.mdp; do
    if grep -qi "^[[:space:]]*define[[:space:]]*=.*POSRES" "$MDP_DIR/$m" \
       && grep -qiE "^[[:space:]]*pcoupl[[:space:]]*=[[:space:]]*(berendsen|c-rescale|parrinello|mttk)" "$MDP_DIR/$m" \
       && ! grep -qi "^[[:space:]]*refcoord.scaling[[:space:]]*=[[:space:]]*com" "$MDP_DIR/$m"; then
        echo "  [FATAL] $m has -DPOSRES + pressure coupling but no"
        echo "          'refcoord-scaling = com'. The barostat would scale atoms"
        echo "          while restraint references stay fixed, and the run would"
        echo "          explode a few ps in."
        echo "          Fix:  echo 'refcoord-scaling = com' >> $MDP_DIR/$m"
        exit 14
    fi
done
echo "  [ok] refcoord-scaling = com present in restrained NPT stages"
echo ""

# ============================ WORKSPACE ============================
mkdir -p "$WORK"
cd "$WORK"
cp -f "$INPUT_PDB" input.pdb
cp -f "$MDP_DIR"/*.mdp .

echo "================================================================"
echo "  SETUP: $SYSNAME"
echo "  force field : $FORCE_FIELD + $WATER_MODEL"
echo "  box         : dodecahedron, $BOX_PADDING nm padding"
echo "  salt        : $SALT_CONC M NaCl"
echo "  accel       : $ACCEL"
echo "================================================================"

# ---- 1. topology ----
if [ -s processed.gro ]; then
    echo "[1/8] pdb2gmx — done, skipping"
else
    echo "[1/8] pdb2gmx — topology (3 chains, NO merge)"
    gmx pdb2gmx -f input.pdb -o processed.gro -p topol.top \
        -ignh -merge no -ff "$FORCE_FIELD" -water "$WATER_MODEL" \
        2>&1 | tee pdb2gmx.log
    [ -s processed.gro ] || { echo "[ERROR] pdb2gmx failed — see pdb2gmx.log"; exit 1; }
fi
echo "      chain topologies: $(ls topol_Protein_chain_*.itp 2>/dev/null | wc -l)"

# ---- 2. box ----
if [ -s boxed.gro ]; then
    echo "[2/8] editconf — done, skipping"
else
    echo "[2/8] editconf — box"
    gmx editconf -f processed.gro -o boxed.gro -c -d "$BOX_PADDING" \
        -bt dodecahedron 2>&1 | tee editconf.log
    [ -s boxed.gro ] || { echo "[ERROR] editconf failed"; exit 2; }
fi

# ---- 3. solvate ----
if [ -s solvated.gro ]; then
    echo "[3/8] solvate — done, skipping"
else
    echo "[3/8] solvate — TIP3P water"
    gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top \
        2>&1 | tee solvate.log
    [ -s solvated.gro ] || { echo "[ERROR] solvate failed"; exit 3; }
fi

# ---- 4. ions ----
if [ -s ionised.gro ]; then
    echo "[4/8] genion — done, skipping"
else
    echo "[4/8] genion — $SALT_CONC M NaCl + neutralise"
    gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 \
        2>&1 | tee grompp_ions.log
    echo "SOL" | gmx genion -s ions.tpr -o ionised.gro -p topol.top \
        -pname NA -nname CL -neutral -conc "$SALT_CONC" 2>&1 | tee genion.log
    [ -s ionised.gro ] || { echo "[ERROR] genion failed"; exit 4; }
fi

# ---- 5. energy minimisation ----
if [ -s em.gro ]; then
    echo "[5/8] EM — done, skipping"
else
    echo "[5/8] energy minimisation (steepest descent)"
    gmx grompp -f em.mdp -c ionised.gro -p topol.top -o em.tpr -maxwarn 2 \
        2>&1 | tee grompp_em.log
    gmx mdrun -deffnm em -ntmpi 1 -ntomp "$NTOMP" 2>&1 | tee mdrun_em.log
    [ -s em.gro ] || { echo "[ERROR] EM failed — see mdrun_em.log"; exit 5; }
fi
grep -E "Potential Energy|Maximum force" em.log | tail -2 | sed 's/^/      /'
echo "      (Fmax below ~5000 kJ/mol/nm is comfortable; Epot=inf means the"
echo "       INPUT has overlapping atoms — check it, do not minimise harder)"

# ---- 6. NVT ----
if ls nvt*.gro >/dev/null 2>&1; then
    echo "[6/8] NVT — done, skipping"
else
    echo "[6/8] NVT (100 ps, restrained)"
    [ -s nvt.tpr ] || gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top \
        -o nvt.tpr -maxwarn 2 2>&1 | tee grompp_nvt.log
    gmx mdrun -deffnm nvt -cpi nvt.cpt $ACCEL -ntmpi 1 -ntomp "$NTOMP" \
        -maxh "$MAXH" 2>&1 | tee -a mdrun_nvt.log
    ls nvt*.gro >/dev/null 2>&1 || {
        echo "[STOPPED] NVT incomplete (walltime?). Resubmit to resume."; exit 6; }
fi

# ---- 7. NPT warm-up ----
if ls npt_warmup*.gro >/dev/null 2>&1; then
    echo "[7/8] NPT warm-up — done, skipping"
else
    echo "[7/8] NPT warm-up (50 ps @ 0.5 fs, restrained)"
    [ -s npt_warmup.tpr ] || gmx grompp -f npt_warmup.mdp -c nvt.gro -r nvt.gro \
        -t nvt.cpt -p topol.top -o npt_warmup.tpr -maxwarn 2 \
        2>&1 | tee grompp_warmup.log
    gmx mdrun -deffnm npt_warmup -cpi npt_warmup.cpt $ACCEL \
        -ntmpi 1 -ntomp "$NTOMP" -maxh "$MAXH" 2>&1 | tee -a mdrun_warmup.log
    ls npt_warmup*.gro >/dev/null 2>&1 || {
        echo "[STOPPED] warm-up incomplete. Resubmit to resume."; exit 7; }
fi

# ---- 8. NPT main ----
if ls npt.gro npt.part*.gro >/dev/null 2>&1; then
    echo "[8/8] NPT — done, skipping"
else
    echo "[8/8] NPT (1 ns, restrained)"
    [ -s npt.tpr ] || gmx grompp -f npt.mdp -c npt_warmup.gro -r npt_warmup.gro \
        -t npt_warmup.cpt -p topol.top -o npt.tpr -maxwarn 2 \
        2>&1 | tee grompp_npt.log
    gmx mdrun -deffnm npt -cpi npt.cpt $ACCEL -ntmpi 1 -ntomp "$NTOMP" \
        -maxh "$MAXH" 2>&1 | tee -a mdrun_npt.log
    ls npt.gro npt.part*.gro >/dev/null 2>&1 || {
        echo "[STOPPED] NPT incomplete. Resubmit the SAME job to resume."; exit 8; }
fi

echo "================================================================"
echo "  SETUP COMPLETE — $SYSNAME"
echo "  Ready for production: qsub job_prod.pbs"
echo "================================================================"
