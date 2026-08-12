#!/usr/bin/env python3
"""
01_generate_jobs.py — write one PBS job per (ligand x target) pair.

Generalised from the project's working HADDOCK 2.5 batch scripts. The bug fixes
that made those scripts reliable are preserved here and marked in the generated
job files, because each one cost a full queue cycle to find:

  1. Result collection does NOT inherit `set -e`. Each copy is independent, so a
     missing file.list cannot abort collection of the structures themselves.
     This was the bug that left work_pbs/ full and results/ empty.
  2. The PBS log directory is created before submission, not by the job.
  3. `set -e` is scoped to setup and docking only, never to post-processing.
  4. The first run_haddock.py call returns non-zero BY DESIGN (it only generates
     run.cns), so `set +e` wraps it. Treating that as failure aborts every job.
  5. Idempotent: a pair whose water-refined results are already collected is
     skipped, so resubmitting the whole batch is always safe.
  6. Job names are sanitised and truncated — PBS rejects long or odd names.

Usage:
    python3 3_docking/01_generate_jobs.py config/hadv_c5_hvr7.yaml
    python3 3_docking/01_generate_jobs.py config/... --target HD5   # one target
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from common import (banner, get_ligand_dir, get_targets, load_config,
                    model_id, write_manifest)


def safe_job_name(target: str, ligand: str, limit: int = 15) -> str:
    raw = f"{target}_{ligand}"
    return re.sub(r"[^A-Za-z0-9_-]", "_", raw)[:limit]


HADDOCK25_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -l select=1:ncpus={ncpu}:mem={mem}
#PBS -l walltime={walltime}
#PBS -j oe
#PBS -o {logs_dir}/{ligand}.log

# =====================================================================
# ENGINE: HADDOCK 2.5   (run.param + run.cns interface)
# Docking: {target} <-> {ligand}
#
# Resubmitting this script is safe: it skips work already collected.
# =====================================================================

module load python/3.10.9

HADDOCK_DIR="{haddock_dir}"
BASE_DIR="{base_dir}"
TARGET_PDB="{base_dir}/{target}.pdb"
AIR_FILE="{base_dir}/{target}_air.tbl"
LIGAND_PDB="{ligand_pdb}"
RESULTS_DIR="{base_dir}/results"
WORK_DIR="{base_dir}/work_pbs"
LIGAND="{ligand}"
NCPU={ncpu}

PYTHON="python3"
RESULT_DIR="$RESULTS_DIR/$LIGAND"

# ---- Phase 1: idempotency ----
if [ -f "$RESULT_DIR/it1/file.list" ] && \\
   [ "$(ls "$RESULT_DIR/it1/"*.pdb 2>/dev/null | wc -l)" -ge {n_water} ]; then
    echo "[SKIP] $LIGAND already complete"
    exit 0
fi

# ---- Phase 2: setup (strict: fail fast on bad config) ----
set -euo pipefail
source "$HADDOCK_DIR/haddock_configure.sh"
mkdir -p "$RESULTS_DIR" "$WORK_DIR"

run_dir="$WORK_DIR/run_$LIGAND"
mkdir -p "$run_dir"
cd "$run_dir"

cat > run.param << PARAM
HADDOCK_DIR=$HADDOCK_DIR
N_COMP=2
PDB_FILE1=$LIGAND_PDB
PDB_FILE2=$TARGET_PDB
PROJECT_DIR=./docking_$LIGAND
RUN_NUMBER=1
AMBIG_TBL=$AIR_FILE
STRUCTURES_0={n_it0}
STRUCTURES_1={n_it1}
WATERREFINE={n_water}
PARAM

# ---- Phase 3: setup pass ----
# NOTE: this call exits non-zero by design; it only writes run.cns.
set +e
$PYTHON "$HADDOCK_DIR/haddock/run_haddock.py"
set -e

if [ ! -d "docking_$LIGAND/run1" ]; then
    echo "[ERROR] HADDOCK setup failed for $LIGAND — no run1 directory" >&2
    exit 11
fi
cd "docking_$LIGAND/run1"

# ---- Phase 4: patch run.cns ----
sed -i \\
    -e "s/structures_0=1000/structures_0={n_it0}/g" \\
    -e "s/structures_1=200/structures_1={n_it1}/g" \\
    -e "s/waterrefine=200/waterrefine={n_water}/g" \\
    -e "s/cpunumber_1=[0-9]*/cpunumber_1=$NCPU/g" \\
    -e "s/clust_meth=.*/clust_meth=FCC;/g" \\
    -e "s/clust_cutoff=.*/clust_cutoff={clust_cutoff};/g" \\
    -e "s/clust_size=.*/clust_size={clust_size};/g" \\
    run.cns

# ---- Phase 5: dock ----
echo "[INFO] docking $LIGAND vs {target} started $(date)"
$PYTHON "$HADDOCK_DIR/haddock/run_haddock.py"
echo "[INFO] docking finished $(date)"

# ---- Phase 6: collect (NON-strict on purpose) ----
# Each copy is independent. A missing file.list must not prevent the
# structures themselves from being collected.
set +e
mkdir -p "$RESULT_DIR/it0" "$RESULT_DIR/it1"

if [ -d structures/it0 ]; then
    cp structures/it0/*.pdb "$RESULT_DIR/it0/" 2>/dev/null
    cp structures/it0/file.list "$RESULT_DIR/it0/" 2>/dev/null
fi

n_water=0
if [ -d structures/it1/water ]; then
    n_water=$(ls structures/it1/water/*.pdb 2>/dev/null | wc -l)
    cp structures/it1/water/*.pdb "$RESULT_DIR/it1/" 2>/dev/null
    cp structures/it1/water/file.list "$RESULT_DIR/it1/" 2>/dev/null
fi
cp -r structures/it1/analysis "$RESULT_DIR/" 2>/dev/null

# ---- Phase 7: validate ----
if [ "$n_water" -lt {n_water} ]; then
    echo "[ERROR] expected {n_water} water-refined PDBs, got $n_water" >&2
    echo "        structures may still be in $run_dir — run 02_collect_results.py" >&2
    exit 12
fi
echo "[SUCCESS] $LIGAND: $n_water water-refined structures"
exit 0
"""


HADDOCK3_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -l select=1:ncpus={ncpu}:mem={mem}
#PBS -l walltime={walltime}
#PBS -j oe
#PBS -o {logs_dir}/{ligand}.log

# =====================================================================
# ENGINE: HADDOCK3   (modular TOML workflow)
# Docking: {target} <-> {ligand}
#
# NOT interchangeable with the 2.5 jobs: different modules, different
# output layout, and scores that should not be pooled with 2.5 scores.
# =====================================================================

BASE_DIR="{base_dir}"
RUN_DIR="$BASE_DIR/work_pbs/run_{ligand}"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

if [ -d "$BASE_DIR/results/{ligand}/haddock3" ]; then
    echo "[SKIP] {ligand} already has HADDOCK3 results"
    exit 0
fi

cat > run.toml << 'TOML'
run_dir = "haddock3_out"
molecules = [
    "{ligand_pdb}",
    "{base_dir}/{target}.pdb",
]

[topoaa]

[rigidbody]
ambig_fname = "{base_dir}/{target}_air.tbl"
sampling = {n_it0}

[caprieval]

[flexref]
ambig_fname = "{base_dir}/{target}_air.tbl"

[emref]
ambig_fname = "{base_dir}/{target}_air.tbl"

[clustfcc]
clust_cutoff = {clust_cutoff}
min_population = {clust_size}

[seletopclusts]

[caprieval]
TOML

echo "[INFO] HADDOCK3 docking {ligand} vs {target} started $(date)"
haddock3 run.toml
rc=$?
echo "[INFO] finished $(date) (exit $rc)"

mkdir -p "$BASE_DIR/results/{ligand}"
cp -r haddock3_out "$BASE_DIR/results/{ligand}/haddock3" 2>/dev/null
exit $rc
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--target", default=None,
                    help="generate jobs for one target only")
    ap.add_argument("--out", default=None,
                    help="base directory (default: docking.base_dir or data/docking)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    dock = cfg["docking"]
    engine = dock["engine"]
    targets = get_targets(cfg)
    lig_dir = get_ligand_dir(cfg)

    if args.target:
        if args.target not in targets:
            sys.exit(f"[FATAL] target '{args.target}' not in config "
                     f"(have: {', '.join(targets)})")
        targets = {args.target: targets[args.target]}

    ligands = sorted(lig_dir.glob("*.pdb"))
    if not ligands:
        sys.exit(f"[FATAL] no ligand PDBs in {lig_dir}")

    base_root = Path(args.out or dock.get("base_dir", "data/docking")).expanduser()
    banner(f"Generating {engine} jobs — "
           f"{len(ligands)} ligand(s) x {len(targets)} target(s)")

    total = 0
    for tname, spec in targets.items():
        base = base_root / tname
        logs = base / "logs"
        for d in (logs, base / "results", base / "work_pbs", base / "jobs"):
            d.mkdir(parents=True, exist_ok=True)

        # stage target PDB + AIR so each job is self-contained
        (base / f"{tname}.pdb").write_bytes(spec["pdb"].read_bytes())
        (base / f"{tname}_air.tbl").write_bytes(spec["air"].read_bytes())

        tmpl = HADDOCK25_TEMPLATE if engine == "haddock2.5" else HADDOCK3_TEMPLATE
        for lig in ligands:
            lid = model_id(lig)
            job = tmpl.format(
                job_name=safe_job_name(tname, lid),
                ncpu=dock["ncpu"], mem=dock["mem"], walltime=dock["walltime"],
                logs_dir=logs, base_dir=base, target=tname,
                ligand=lid, ligand_pdb=lig.resolve(),
                haddock_dir=Path(dock.get("haddock_dir",
                                          "~/software/haddock2.5-2025-08")).expanduser(),
                n_it0=dock["n_it0"], n_it1=dock["n_it1"], n_water=dock["n_water"],
                clust_cutoff=dock["clust_cutoff"], clust_size=dock["clust_size"],
            )
            dest = base / "jobs" / f"job_{lid}.pbs"
            dest.write_text(job)
            dest.chmod(0o755)
            total += 1

        print(f"  [{tname}] {len(ligands)} jobs -> {base/'jobs'}")

    print()
    print(f"  {total} job(s) written, engine = {engine}")
    print()
    print("  Submit:")
    for tname in targets:
        print(f"    for f in {base_root/tname/'jobs'}/*.pbs; do qsub \"$f\"; done")
    print()
    print("  If results/ ends up empty but work_pbs/ has structures:")
    print("    python3 3_docking/02_collect_results.py <config>")

    write_manifest(base_root, "generate_jobs", args.config, ligands,
                   extra={"engine": engine, "targets": list(targets),
                          "n_jobs": total, "sampling": {
                              "it0": dock["n_it0"], "it1": dock["n_it1"],
                              "water": dock["n_water"]}})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
