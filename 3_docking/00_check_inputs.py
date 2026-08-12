#!/usr/bin/env python3
"""
00_check_inputs.py — verify everything docking needs, before queuing anything.

Docking runs are long and are submitted in bulk, so a missing AIR file or an
unset CNS path costs a whole queue cycle to discover. Every check here maps to
a failure that actually happened while building this pipeline.

Checks
  1. HADDOCK install present, and which engine it is
  2. CNS executable present (HADDOCK cannot run without it)
  3. Every target has both a PDB and an AIR restraint file
  4. AIR files parse and contain restraints
  5. Ligand structures exist and are SINGLE CHAIN
  6. No filename characters that break PBS job names

Usage:
    python3 3_docking/00_check_inputs.py config/hadv_c5_hvr7.yaml
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

from common import banner, get_ligand_dir, get_targets, load_config


def check_haddock(cfg) -> tuple[bool, str]:
    dock = cfg["docking"]
    engine = dock["engine"]
    ok = True

    if engine == "haddock2.5":
        hdir = Path(dock.get("haddock_dir", "~/software/haddock2.5-2025-08")).expanduser()
        if not hdir.is_dir():
            print(f"  [FAIL] HADDOCK 2.5 not found at {hdir}")
            print("         Set docking.haddock_dir in the config.")
            print("         HADDOCK 2.5 is requested from the Bonvin lab; licences")
            print("         are personal and cannot be redistributed.")
            return False, engine
        cfgsh = hdir / "haddock_configure.sh"
        if not cfgsh.is_file():
            print(f"  [FAIL] {cfgsh} missing — install looks incomplete")
            ok = False
        else:
            print(f"  [ok] HADDOCK 2.5 at {hdir}")
        run_py = hdir / "haddock" / "run_haddock.py"
        if not run_py.is_file():
            print(f"  [FAIL] {run_py} missing")
            ok = False
    else:
        from shutil import which
        exe = which("haddock3")
        if not exe:
            print("  [FAIL] haddock3 not on PATH")
            print("         pip install haddock3   (https://github.com/haddocking/haddock3)")
            return False, engine
        print(f"  [ok] HADDOCK3 at {exe}")

    return ok, engine


def check_cns(cfg) -> bool:
    """
    CNS is the computational engine underneath both HADDOCK versions. It is free
    for non-profit use but must be registered for separately at cns-online.org —
    it is not bundled with HADDOCK and not redistributable here.
    """
    dock = cfg["docking"]
    cns = dock.get("cns_exe")
    if cns:
        p = Path(cns).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            print(f"  [ok] CNS executable {p}")
            return True
        print(f"  [FAIL] docking.cns_exe set but not executable: {p}")
        return False

    env = os.environ.get("CNS_EXE") or os.environ.get("CNS_SOLVE")
    if env and Path(env).exists():
        print(f"  [ok] CNS from environment: {env}")
        return True

    print("  [WARN] CNS not verified. HADDOCK will fail at run time without it.")
    print("         Set docking.cns_exe in the config, or source cns_solve_env.")
    print("         CNS is free for non-profit use: http://cns-online.org/v1.3/")
    return True   # warn, do not block — it may be resolved by haddock_configure.sh


def count_restraints(path: Path) -> int:
    try:
        text = path.read_text()
    except OSError:
        return 0
    return len(re.findall(r"\bassign\b", text, flags=re.IGNORECASE))


def check_targets(cfg) -> bool:
    targets = get_targets(cfg)
    ok = True
    print(f"  panel: {len(targets)} target(s)")
    for name, spec in targets.items():
        pdb, air = spec["pdb"], spec["air"]
        if not pdb.is_file():
            print(f"  [FAIL] {name}: PDB not found — {pdb}")
            ok = False
            continue
        if not air.is_file():
            print(f"  [FAIL] {name}: AIR file not found — {air}")
            print("         AIR restraints are an INPUT you provide. Generate them")
            print("         from your interface residues (e.g. the Bonvin lab web")
            print("         tool) — they encode which residues should contact.")
            ok = False
            continue
        n = count_restraints(air)
        if n == 0:
            print(f"  [FAIL] {name}: {air.name} contains no 'assign' restraints")
            ok = False
        else:
            chains = {l[21] for l in pdb.read_text().splitlines()
                      if l.startswith("ATOM")}
            print(f"  [ok] {name}: {pdb.name} (chains {''.join(sorted(chains))}), "
                  f"{air.name} ({n} restraints)")
    return ok


def check_ligands(cfg) -> bool:
    lig_dir = get_ligand_dir(cfg)
    if not lig_dir.is_dir():
        print(f"  [FAIL] ligand directory not found: {lig_dir}")
        print("         Run Tool 1 stage 07 to produce single-chain structures,")
        print("         or set paths.models_docking in the config.")
        return False

    pdbs = sorted(lig_dir.glob("*.pdb"))
    if not pdbs:
        print(f"  [FAIL] no .pdb files in {lig_dir}")
        return False

    ok = True
    multi = []
    for p in pdbs:
        chains = {l[21] for l in p.read_text().splitlines() if l.startswith("ATOM")}
        if len(chains) > 1:
            multi.append((p.name, "".join(sorted(chains))))
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", p.name):
            print(f"  [FAIL] {p.name}: filename has characters that break PBS job names")
            ok = False

    print(f"  [ok] {len(pdbs)} ligand structure(s) in {lig_dir}")
    if multi:
        ok = False
        print(f"  [FAIL] {len(multi)} structure(s) have more than one chain:")
        for name, ch in multi[:5]:
            print(f"         {name} (chains {ch})")
        print("         HADDOCK treats each input as ONE rigid body. A three-chain")
        print("         trimer would be docked as three separate molecules.")
        print("         Use the single-chain output from Tool 1 stage 07.")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    args = ap.parse_args()

    cfg = load_config(args.config)
    banner("Docking pre-flight")

    print("\n[1] HADDOCK")
    hd_ok, engine = check_haddock(cfg)

    print("\n[2] CNS")
    cns_ok = check_cns(cfg)

    print("\n[3] Targets")
    tg_ok = check_targets(cfg)

    print("\n[4] Ligand structures")
    lg_ok = check_ligands(cfg)

    dock = cfg["docking"]
    print("\n[5] Sampling")
    print(f"  engine   : {engine}")
    print(f"  it0      : {dock['n_it0']} rigid-body models")
    print(f"  it1      : {dock['n_it1']} semi-flexible")
    print(f"  water    : {dock['n_water']} final water-refined")
    print(f"  resources: {dock['ncpu']} cpus, {dock['mem']}, {dock['walltime']}")

    all_ok = hd_ok and cns_ok and tg_ok and lg_ok
    print()
    banner("READY — generate jobs with 01_generate_jobs.py" if all_ok
           else "NOT READY — fix the [FAIL] items above")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
