#!/usr/bin/env python3
"""
03_region_sasa.py — solvent-accessible surface area per region.

Why this matters for humanisation: an engineered region only evades an antibody
if it is still presented on the capsid surface. A graft that folds inward, or
one that balloons outward relative to wild type, has changed its antigenic
presentation even when pLDDT and geometry look fine. SASA is the cheapest way
to see that before committing to docking.

This is a comparative measure, not an absolute one. The number to look at is
the ratio against the wild-type reference, which is why a reference model can
be named in the config (validation.reference) or passed with --reference.

Adapted from the project's original HVR SASA script, with three changes:
  * region boundaries come from the config, not a hard-coded dict
  * residue numbering is folded onto protomer 1, so it works on both
    three-chain models and merged single-chain models
  * output is per-chain as well as summed, so an asymmetric trimer is visible

Usage:
    python3 03_region_sasa.py config/hadv_c5_hvr7.yaml [--models DIR]
                              [--reference wild_type]

Writes: reports/validation/sasa.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley
except ImportError:
    sys.exit("Biopython missing.  pip install -r requirements.txt")

from common import (banner, get_regions, load_config, model_id,
                    normalise_resnum, protomer_length, write_manifest)


def load_structure(path: Path):
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    return parser.get_structure("m", str(path))


def region_sasa(structure, bounds, plen) -> tuple[float, dict]:
    """Total SASA in a region, plus a per-chain breakdown."""
    lo, hi = bounds
    total = 0.0
    per_chain: dict[str, float] = {}
    for model in structure:
        for chain in model:
            acc = 0.0
            for res in chain:
                if res.id[0] != " ":         # skip waters/hetero
                    continue
                if lo <= normalise_resnum(res.id[1], plen) <= hi:
                    acc += getattr(res, "sasa", 0.0)
            if acc:
                per_chain[chain.id] = per_chain.get(chain.id, 0.0) + acc
            total += acc
        break                                 # first model only
    return total, per_chain


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--models", default=None)
    ap.add_argument("--reference", default=None,
                    help="model id to use as the comparison baseline "
                         "(default: validation.reference in config)")
    ap.add_argument("--n-points", type=int, default=100,
                    help="Shrake-Rupley sphere points; 960 is slower but smoother")
    ap.add_argument("--out", default="reports/validation")
    args = ap.parse_args()

    cfg = load_config(args.config)
    regions = get_regions(cfg)
    plen = protomer_length(cfg)
    reference = args.reference or (cfg.get("validation") or {}).get("reference")

    model_dir = Path(args.models or (cfg.get("paths") or {}).get("models", "data/models"))
    if not model_dir.is_dir():
        sys.exit(f"[FATAL] model directory not found: {model_dir}")

    files = sorted([p for p in model_dir.iterdir() if p.suffix in (".pdb", ".cif")])
    if not files:
        sys.exit(f"[FATAL] no models in {model_dir}")

    banner(f"Region SASA — {len(files)} model(s)")
    print(f"  regions: {', '.join(regions)}")
    if reference:
        print(f"  reference: {reference}")
    print("  (Shrake-Rupley; this is the slow step — a few seconds per model)")
    print()

    sr = ShrakeRupley(probe_radius=1.40, n_points=args.n_points)
    rows = []

    for path in files:
        mid = model_id(path)
        try:
            structure = load_structure(path)
            sr.compute(structure, level="R")
        except Exception as exc:
            print(f"  [ERROR] {mid}: {exc}")
            rows.append({"model": mid, "error": str(exc)[:120]})
            continue

        row = {"model": mid}
        combined = 0.0
        for name, bounds in regions.items():
            total, per_chain = region_sasa(structure, bounds, plen)
            row[f"sasa_{name}"] = round(total, 2)
            combined += total
            if per_chain:
                spread = (max(per_chain.values()) - min(per_chain.values()))
                row[f"sasa_{name}_chain_spread"] = round(spread, 2)
        row["sasa_all_regions"] = round(combined, 2)

        chains = sorted({c.id for m in structure for c in m})
        row["n_chains"] = len(chains)
        rows.append(row)
        print(f"  {mid:<40} total region SASA = {row['sasa_all_regions']:>10.1f} A^2")

    # ratios against the reference, if we have one
    ref_row = next((r for r in rows if r.get("model") == reference), None)
    if ref_row:
        for row in rows:
            for name in regions:
                key, ref_val = f"sasa_{name}", ref_row.get(f"sasa_{name}")
                if isinstance(ref_val, (int, float)) and ref_val > 0 \
                        and isinstance(row.get(key), (int, float)):
                    row[f"ratio_{name}"] = round(row[key] / ref_val, 3)
    elif reference:
        print(f"\n  [note] reference '{reference}' not found among models — "
              f"ratios skipped")

    fieldnames, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "sasa.csv"
    with open(dest, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"  wrote {dest}")
    if ref_row:
        print(f"  ratio_* columns are relative to {reference}; "
              f"~1.0 means presentation is preserved")
    write_manifest(out_dir, "sasa", args.config, files,
                   extra={"reference": reference, "n_points": args.n_points})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
