#!/usr/bin/env python3
"""
02_collect_results.py — rescue results from work_pbs into results/.

Why this exists as a separate stage: in the original run, HADDOCK completed
successfully but `results/` stayed empty while `work_pbs/` held every structure.
The collection block had inherited `set -e`, so the first failed copy silently
aborted the rest. Jobs reported failure despite the science being finished.

The generated jobs now collect their own output, but this stage stays because:
  * it recovers older runs that predate the fix
  * it recovers jobs killed after docking but before collection
  * it is idempotent and safe to run at any time

Usage:
    python3 3_docking/02_collect_results.py config/hadv_c5_hvr7.yaml
    python3 3_docking/02_collect_results.py config/... --target HD5
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from common import banner, get_targets, load_config, write_manifest


def collect_25(work_dir: Path, results_dir: Path) -> list[dict]:
    """HADDOCK 2.5 layout: work_pbs/run_<lig>/docking_<lig>/run1/structures/"""
    rows = []
    for run_d in sorted(work_dir.glob("run_*")):
        lig = run_d.name[len("run_"):]
        inner = run_d / f"docking_{lig}" / "run1" / "structures"
        if not inner.is_dir():
            # tolerate a differently-named project dir
            cand = list(run_d.glob("docking_*/run1/structures"))
            if not cand:
                rows.append({"ligand": lig, "it0": 0, "water": 0,
                             "note": "no structures dir"})
                continue
            inner = cand[0]

        out = results_dir / lig
        (out / "it0").mkdir(parents=True, exist_ok=True)
        (out / "it1").mkdir(parents=True, exist_ok=True)

        n_it0 = n_water = 0

        src0 = inner / "it0"
        if src0.is_dir():
            for p in src0.glob("*.pdb"):
                shutil.copy2(p, out / "it0" / p.name)
                n_it0 += 1
            fl = src0 / "file.list"
            if fl.is_file():
                shutil.copy2(fl, out / "it0" / "file.list")

        srcw = inner / "it1" / "water"
        if srcw.is_dir():
            for p in srcw.glob("*.pdb"):
                shutil.copy2(p, out / "it1" / p.name)
                n_water += 1
            fl = srcw / "file.list"
            if fl.is_file():
                shutil.copy2(fl, out / "it1" / "file.list")

        # cluster analysis output — needed for cluster-level ranking
        ana = inner / "it1" / "analysis"
        if ana.is_dir():
            dest = out / "analysis"
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(ana, dest)

        has_list = (out / "it1" / "file.list").is_file()
        rows.append({"ligand": lig, "it0": n_it0, "water": n_water,
                     "file_list": has_list,
                     "note": "" if has_list else "no it1/file.list"})
    return rows


def collect_3(work_dir: Path, results_dir: Path) -> list[dict]:
    """HADDOCK3 layout: work_pbs/run_<lig>/haddock3_out/<NN_module>/"""
    rows = []
    for run_d in sorted(work_dir.glob("run_*")):
        lig = run_d.name[len("run_"):]
        out_dir = run_d / "haddock3_out"
        if not out_dir.is_dir():
            rows.append({"ligand": lig, "modules": 0, "note": "no haddock3_out"})
            continue
        dest = results_dir / lig / "haddock3"
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(out_dir, dest)
        mods = [p.name for p in sorted(dest.iterdir()) if p.is_dir()]
        rows.append({"ligand": lig, "modules": len(mods),
                     "note": ",".join(mods[:4])})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--target", default=None)
    ap.add_argument("--base", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    dock = cfg["docking"]
    engine = dock["engine"]
    targets = get_targets(cfg)
    if args.target:
        if args.target not in targets:
            sys.exit(f"[FATAL] unknown target '{args.target}'")
        targets = {args.target: targets[args.target]}

    base_root = Path(args.base or dock.get("base_dir", "data/docking")).expanduser()
    banner(f"Collecting {engine} results")

    all_rows = []
    for tname in targets:
        base = base_root / tname
        work = base / "work_pbs"
        res = base / "results"
        if not work.is_dir():
            print(f"  [{tname}] no work_pbs — nothing to collect")
            continue
        res.mkdir(parents=True, exist_ok=True)

        rows = (collect_25 if engine == "haddock2.5" else collect_3)(work, res)
        print(f"\n  [{tname}]")
        for r in rows:
            if engine == "haddock2.5":
                flag = "ok " if r.get("water", 0) > 0 else "EMPTY"
                print(f"    {flag} {r['ligand']:<40} it0={r.get('it0',0):<5} "
                      f"water={r.get('water',0):<5} {r.get('note','')}")
            else:
                print(f"    {r['ligand']:<40} modules={r.get('modules',0)} "
                      f"{r.get('note','')}")
            r["target"] = tname
        all_rows += rows

    ok = sum(1 for r in all_rows
             if (r.get("water", 0) > 0 or r.get("modules", 0) > 0))
    print()
    print(f"  {ok}/{len(all_rows)} pair(s) have collected results")
    if ok < len(all_rows):
        print("  Empty ones: check logs/<ligand>.log — docking may have failed,")
        print("  or the job may still be running.")

    write_manifest(base_root, "collect", args.config, [],
                   extra={"engine": engine, "n_ok": ok, "n_total": len(all_rows)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
