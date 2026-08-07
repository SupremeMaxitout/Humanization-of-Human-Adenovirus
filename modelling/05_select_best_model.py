#!/usr/bin/env python3
"""
05_select_best_model.py — pick the best seed per candidate and record ALL
confidence metrics.

The original kept only the winning model and printed scores to the terminal.
Here every metric is captured to a spreadsheet, because these ARE the
validation data for Tool 2:

  AlphaFold3  (*_summary_confidences_*.json)
      ranking_score, ptm, iptm, fraction_disordered, has_clash, num_recycles
  ColabFold   (*_scores_rank_*.json)
      plddt (per-residue -> mean), ptm, iptm, max_pae

Rejects models flagged has_clash when selection.reject_if_has_clash is set.

Usage:
    python3 05_select_best_model.py config/hadv_c5_hvr7.yaml
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, write_manifest, banner, ensure_dirs

ap = argparse.ArgumentParser()
ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
args = ap.parse_args()

cfg = load_config(args.config)
sel = cfg["selection"]
ensure_dirs(cfg, "best_models", "reports")

extracted = Path(cfg["paths"]["extracted"])
best_dir = Path(cfg["paths"]["best_models"])
rank_by = sel.get("rank_by", "ranking_score")

banner(f"Selecting best model per candidate (rank_by = {rank_by})")

rows, copied, skipped = [], 0, []

for d in sorted(p for p in extracted.iterdir() if p.is_dir()):
    af3_jsons = sorted(d.glob("*summary_confidences*.json"))
    cf_jsons = sorted(d.glob("*scores_rank*.json"))

    cands = []

    # ---- AlphaFold3 ----
    for j in af3_jsons:
        try:
            data = json.loads(j.read_text())
        except Exception as e:
            print(f"  [warn] unreadable {j.name}: {e}")
            continue
        idx = j.stem.split("_")[-1]
        cands.append({
            "seed": idx,
            "engine": "alphafold3",
            "score": data.get(rank_by, data.get("ranking_score",
                        data.get("ptm", data.get("iptm")))),
            "ranking_score": data.get("ranking_score"),
            "ptm": data.get("ptm"),
            "iptm": data.get("iptm"),
            "fraction_disordered": data.get("fraction_disordered"),
            "has_clash": data.get("has_clash"),
            "model_glob": f"*model_{idx}.cif",
        })

    # ---- ColabFold ----
    for j in cf_jsons:
        try:
            data = json.loads(j.read_text())
        except Exception:
            continue
        plddt = data.get("plddt") or []
        mean_plddt = round(sum(plddt) / len(plddt), 2) if plddt else None
        rank_tag = j.stem.split("rank_")[-1].split("_")[0] if "rank_" in j.stem else "001"
        cands.append({
            "seed": rank_tag,
            "engine": "colabfold",
            "score": data.get("ptm") if rank_by in ("ranking_score", "ptm")
                     else data.get("iptm"),
            "ranking_score": None,
            "ptm": data.get("ptm"),
            "iptm": data.get("iptm"),
            "mean_plddt": mean_plddt,
            "max_pae": data.get("max_pae"),
            "has_clash": None,
            "model_glob": f"*rank_{rank_tag}*.pdb",
        })

    if not cands:
        skipped.append(f"{d.name} (no confidence JSON)")
        continue

    usable = [c for c in cands if c["score"] is not None]
    if sel.get("reject_if_has_clash", True):
        clean = [c for c in usable if not c.get("has_clash")]
        if clean:
            usable = clean
        elif usable:
            print(f"  [warn] {d.name}: every seed flagged has_clash")

    if not usable:
        skipped.append(f"{d.name} (no usable scores)")
        continue

    best = max(usable, key=lambda c: c["score"])

    if best["score"] < sel.get("min_ranking_score", 0.0):
        skipped.append(f"{d.name} (score {best['score']:.3f} below floor)")
        continue

    hits = sorted(d.glob(best["model_glob"]))
    if not hits:
        skipped.append(f"{d.name} (model file {best['model_glob']} missing)")
        continue

    src = hits[0]
    dst = best_dir / f"{d.name}{src.suffix}"
    shutil.copy2(src, dst)
    copied += 1

    rows.append({
        "candidate": d.name,
        "engine": best["engine"],
        "best_seed": best["seed"],
        "selected_by": rank_by,
        "score": best["score"],
        "ranking_score": best.get("ranking_score"),
        "ptm": best.get("ptm"),
        "iptm": best.get("iptm"),
        "mean_plddt": best.get("mean_plddt"),
        "fraction_disordered": best.get("fraction_disordered"),
        "has_clash": best.get("has_clash"),
        "n_seeds": len(cands),
        "source_file": src.name,
        "output_file": dst.name,
    })
    print(f"  [ok] {d.name}: seed {best['seed']}, {rank_by}={best['score']:.4f}")

df = pd.DataFrame(rows)
out = Path(cfg["paths"]["reports"]) / "model_selection.xlsx"
if not df.empty:
    df = df.sort_values("score", ascending=False)
    df.to_excel(out, index=False)

print(f"\n  selected : {copied}")
print(f"  skipped  : {len(skipped)}")
for s in skipped:
    print(f"    - {s}")
if not df.empty:
    print(f"\n  models -> {best_dir}")
    print(f"  metrics -> {out}")
    print(f"  score range: {df['score'].min():.3f} - {df['score'].max():.3f}")

write_manifest(cfg, "05_select_best_model",
               inputs={"extracted_dir": str(extracted)},
               outputs={"best_models": str(best_dir), "metrics": str(out),
                        "n_selected": copied},
               extra={"rank_by": rank_by, "skipped": skipped})
