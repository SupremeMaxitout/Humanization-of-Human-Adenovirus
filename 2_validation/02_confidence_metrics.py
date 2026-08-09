#!/usr/bin/env python3
"""
02_confidence_metrics.py — predictor confidence, global and per-region.

The point of this stage: a model can have an excellent global score while the
grafted region specifically is a low-confidence, floppy mess. Since the graft is
the entire subject of the study, region-restricted pLDDT matters more than the
overall number, and a mean-only report will hide exactly the failure you care
about.

Metrics collected
  mean_plddt          whole model, from the B-factor column (AF3 and ColabFold
                      both write pLDDT there, 0-100)
  region_plddt        pLDDT restricted to each configured region
  target_region_plddt pLDDT of the region being humanised — the decisive number
  ptm / iptm          fold and interface confidence (iptm = trimer assembly)
  ranking_score       AF3's own composite ranking
  fraction_disordered AF3 estimate of disordered content
  has_clash           AF3's internal clash flag
  interface_pae       mean PAE between chain pairs, where a PAE matrix exists

Supports both AlphaFold3 server output and ColabFold output. AF3 writes
*_summary_confidences_*.json and *_confidences_*.json; ColabFold writes
*_scores_rank_*.json with plddt/pae/ptm/iptm. Missing files are reported as
blanks rather than errors, so a partial set still produces a usable table.

Usage:
    python3 02_confidence_metrics.py config/hadv_c5_hvr7.yaml \
        [--models DIR] [--json-dir DIR]

Writes: reports/validation/confidence.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

try:
    from Bio.PDB import PDBParser, MMCIFParser
except ImportError:
    sys.exit("Biopython missing.  pip install -r requirements.txt")

from common import (banner, get_regions, get_target_region, load_config,
                    model_id, normalise_resnum, protomer_length, write_manifest)


# --------------------------------------------------------------------------
# pLDDT from the structure file
# --------------------------------------------------------------------------

def per_residue_plddt(path: Path) -> dict[tuple[str, int], float]:
    """
    Mean B-factor per residue = pLDDT.

    Works for AF3 mmCIF and ColabFold PDB alike. ColabFold "relaxed" files keep
    pLDDT in the B-factor column too, so this is predictor-agnostic.
    """
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    model = parser.get_structure("m", str(path))[0]
    out = {}
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            vals = [a.get_bfactor() for a in res if a.element != "H"]
            if vals:
                out[(chain.id, res.id[1])] = float(np.mean(vals))
    return out


def region_mean(plddt: dict, bounds: tuple[int, int], plen: int | None) -> float | str:
    lo, hi = bounds
    vals = [v for (_, num), v in plddt.items()
            if lo <= normalise_resnum(num, plen) <= hi]
    return round(float(np.mean(vals)), 2) if vals else ""


# --------------------------------------------------------------------------
# confidence JSONs
# --------------------------------------------------------------------------

def find_json(json_dir: Path, mid: str) -> list[Path]:
    """Candidate confidence JSONs for a model, across both predictors."""
    if not json_dir or not json_dir.is_dir():
        return []
    stem = mid.split("_model")[0]
    hits = []
    for pattern in ("*summary_confidences*.json", "*scores_rank*.json",
                    "*confidences*.json"):
        hits += [p for p in json_dir.rglob(pattern) if stem in p.name or stem in str(p.parent)]
    # de-duplicate, preserve order
    seen, out = set(), []
    for p in hits:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def parse_confidence(paths: list[Path]) -> dict:
    """Pull whatever recognised metrics exist across the given JSONs."""
    got = {}
    pae_matrix = None
    for path in paths:
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception:
            continue
        if not isinstance(data, dict):
            continue

        for key in ("ptm", "iptm", "ranking_score", "fraction_disordered",
                    "has_clash", "num_recycles"):
            if key in data and key not in got:
                val = data[key]
                if isinstance(val, (int, float, bool)):
                    got[key] = val

        # ColabFold puts a flat plddt list here; AF3 uses atom_plddts
        if "pae" in data and pae_matrix is None:
            try:
                pae_matrix = np.asarray(data["pae"], dtype=float)
            except Exception:
                pass
        if "max_pae" in data and "max_pae" not in got:
            got["max_pae"] = data["max_pae"]

        # AF3 chain-pair iptm: mean of off-diagonal entries = interface confidence
        if "chain_pair_iptm" in data and "chain_pair_iptm_mean" not in got:
            try:
                mat = np.asarray(data["chain_pair_iptm"], dtype=float)
                off = mat[~np.eye(mat.shape[0], dtype=bool)]
                if off.size:
                    got["chain_pair_iptm_mean"] = round(float(np.mean(off)), 4)
            except Exception:
                pass

    if pae_matrix is not None and pae_matrix.ndim == 2:
        got["pae_mean"] = round(float(np.mean(pae_matrix)), 3)
        got["pae_max"] = round(float(np.max(pae_matrix)), 3)
    return got


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--models", default=None)
    ap.add_argument("--json-dir", default=None,
                    help="directory containing predictor confidence JSONs "
                         "(searched recursively)")
    ap.add_argument("--out", default="reports/validation")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths_cfg = cfg.get("paths") or {}
    regions = get_regions(cfg)
    target_region = get_target_region(cfg)
    plen = protomer_length(cfg)

    model_dir = Path(args.models or paths_cfg.get("models", "data/models"))
    json_dir = Path(args.json_dir or paths_cfg.get("predictions", "data/predictions"))

    if not model_dir.is_dir():
        sys.exit(f"[FATAL] model directory not found: {model_dir}")

    files = sorted([p for p in model_dir.iterdir() if p.suffix in (".pdb", ".cif")])
    if not files:
        sys.exit(f"[FATAL] no models in {model_dir}")

    banner(f"Confidence metrics — {len(files)} model(s)")
    if plen:
        print(f"  protomer length {plen} (residue numbers folded onto protomer 1)")
    if target_region:
        print(f"  target region: {target_region} {regions.get(target_region)}")
    if not json_dir.is_dir():
        print(f"  [note] no JSON dir at {json_dir} — pLDDT only, no pTM/ipTM/PAE")
    print()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []

    for path in files:
        mid = model_id(path)
        row = {"model": mid}
        try:
            plddt = per_residue_plddt(path)
        except Exception as exc:
            print(f"  [ERROR] {mid}: {exc}")
            rows.append({"model": mid, "error": str(exc)[:120]})
            continue

        vals = list(plddt.values())
        row["mean_plddt"] = round(float(np.mean(vals)), 2) if vals else ""
        row["min_residue_plddt"] = round(float(np.min(vals)), 2) if vals else ""
        row["n_residues"] = len(vals)

        for name, bounds in regions.items():
            row[f"plddt_{name}"] = region_mean(plddt, bounds, plen)

        if target_region and target_region in regions:
            row["target_region"] = target_region
            row["target_region_plddt"] = row.get(f"plddt_{target_region}", "")

        row.update(parse_confidence(find_json(json_dir, mid)))
        rows.append(row)

        tgt = row.get("target_region_plddt", "")
        print(f"  {mid:<40} mean pLDDT={row['mean_plddt']:<7} "
              f"{target_region or 'region'}={tgt}")

    # union of keys, stable order
    fieldnames, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)

    dest = out_dir / "confidence.csv"
    with open(dest, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"  wrote {dest}")
    write_manifest(out_dir, "confidence", args.config, files,
                   extra={"target_region": target_region})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
