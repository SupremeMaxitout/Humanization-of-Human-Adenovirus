#!/usr/bin/env python3
"""
03_score_report.py — cluster-level scoring, ranked against the wild type.

WHY CLUSTERS, NOT THE SINGLE BEST STRUCTURE
-------------------------------------------
A docking run generates hundreds of models. The single lowest-scoring one is the
tail of a distribution and is largely noise: rerun with a different seed and a
different structure wins. What is reproducible is a *cluster* — a group of
models that converged on the same binding mode. HADDOCK's own recommendation is
to rank by the average score of the top members of each cluster, and to treat
cluster size as evidence: a large, well-populated cluster means the sampling
kept finding that mode.

So the primary ranking here is:
    cluster score = mean HADDOCK score of the best N members of cluster 1
with cluster size and the score spread reported alongside, because a cluster of
4 with a wide spread is much weaker evidence than a cluster of 40 with a tight
one.

WHAT THE SCORE IS AND IS NOT
----------------------------
HADDOCK score = 1.0*Evdw + 0.2*Eelec + 1.0*Edesolv + 0.1*Eair
It is an empirical scoring function for ranking poses of ONE complex. It is not
a binding affinity, not a Kd, and not calibrated across different targets. A
score of -120 against HD5 and -120 against Factor X do not mean equal binding.

Comparisons in this report are therefore always made WITHIN one target, as a
delta against the wild type. That is the only comparison the score supports.

Usage:
    python3 3_docking/03_score_report.py config/hadv_c5_hvr7.yaml
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas missing.  pip install -r requirements.txt")

from common import (banner, get_targets, haddock_score, load_config,
                    parse_file_list, parse_haddock_pdb_energies, write_manifest)

TOP_N = 4     # members averaged per cluster — HADDOCK's usual convention


def read_clusters(result_dir: Path) -> dict[int, list[str]]:
    """
    Parse HADDOCK 2.5 cluster membership.

    analysis/cluster.out lines look like:
        Cluster 1 -> 12 45 3 88 ...
    """
    out: dict[int, list[str]] = {}
    for name in ("cluster.out", "clusters.out"):
        path = result_dir / "analysis" / name
        if path.is_file():
            for line in path.read_text().splitlines():
                m = re.match(r"\s*Cluster\s+(\d+)\s*->\s*(.*)", line, re.I)
                if m:
                    out[int(m.group(1))] = m.group(2).split()
            if out:
                return out
    return out


def score_ligand(result_dir: Path) -> dict:
    """Scores for one ligand-target pair, at both structure and cluster level."""
    row: dict = {}
    it1 = result_dir / "it1"
    if not it1.is_dir():
        return {"note": "no it1 directory"}

    ranked = parse_file_list(it1 / "file.list")

    # fall back to reading energies out of the PDB headers
    if not ranked:
        pdbs = sorted(it1.glob("*.pdb"))
        tmp = []
        for p in pdbs:
            sc = haddock_score(parse_haddock_pdb_energies(p))
            if sc is not None:
                tmp.append((p.name, sc))
        ranked = sorted(tmp, key=lambda x: x[1])

    if not ranked:
        return {"note": "no scores found"}

    scores = np.array([s for _, s in ranked], dtype=float)
    row["n_structures"] = len(scores)
    row["best_score"] = round(float(scores.min()), 3)
    row["top4_mean"] = round(float(np.mean(np.sort(scores)[:TOP_N])), 3)
    row["top4_sd"] = round(float(np.std(np.sort(scores)[:TOP_N])), 3)

    # energy components, averaged over the best TOP_N structures
    comps: dict[str, list[float]] = {}
    for name, _ in ranked[:TOP_N]:
        vals = parse_haddock_pdb_energies(it1 / name)
        for k, v in vals.items():
            comps.setdefault(k, []).append(v)
    for k, vs in comps.items():
        row[k] = round(float(np.mean(vs)), 2)

    clusters = read_clusters(result_dir)
    if clusters:
        by_name = dict(ranked)
        best = None
        for cid, members in clusters.items():
            member_scores = []
            for m in members:
                for fname, sc in by_name.items():
                    if re.search(rf"_{re.escape(m)}w?\.pdb$", fname):
                        member_scores.append(sc)
                        break
            if not member_scores:
                continue
            member_scores.sort()
            mean_top = float(np.mean(member_scores[:TOP_N]))
            if best is None or mean_top < best["cluster_score"]:
                best = {"cluster_id": cid,
                        "cluster_size": len(members),
                        "cluster_score": round(mean_top, 3),
                        "cluster_sd": round(float(np.std(member_scores[:TOP_N])), 3)}
        if best:
            row.update(best)
        row["n_clusters"] = len(clusters)
    else:
        row["n_clusters"] = 0
        row["note"] = "no cluster analysis — ranking on top4_mean instead"

    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", default="reports/docking")
    args = ap.parse_args()

    cfg = load_config(args.config)
    dock = cfg["docking"]
    engine = dock["engine"]
    targets = get_targets(cfg)
    reference = (cfg.get("validation") or {}).get("reference", "wild_type")
    base_root = Path(args.base or dock.get("base_dir", "data/docking")).expanduser()

    banner(f"Scoring — engine {engine}, reference {reference}")

    rows = []
    for tname in targets:
        res_dir = base_root / tname / "results"
        if not res_dir.is_dir():
            print(f"  [{tname}] no results directory — skipping")
            continue
        for lig_dir in sorted(p for p in res_dir.iterdir() if p.is_dir()):
            row = {"target": tname, "ligand": lig_dir.name, "engine": engine}
            row.update(score_ligand(lig_dir))
            rows.append(row)

    if not rows:
        sys.exit("[FATAL] no results found — run 02_collect_results.py first")

    df = pd.DataFrame(rows)

    # primary metric: cluster score where available, else top4 mean
    df["rank_metric"] = df.get("cluster_score", pd.Series(dtype=float))
    if "top4_mean" in df.columns:
        df["rank_metric"] = df["rank_metric"].fillna(df["top4_mean"])

    # delta vs reference, computed WITHIN each target only
    df["delta_vs_reference"] = np.nan
    for tname, grp in df.groupby("target"):
        ref = grp.loc[grp["ligand"] == reference, "rank_metric"]
        if not ref.empty and pd.notna(ref.iloc[0]):
            df.loc[grp.index, "delta_vs_reference"] = (
                grp["rank_metric"] - float(ref.iloc[0])).round(3)
        else:
            print(f"  [note] {tname}: reference '{reference}' has no score — "
                  f"deltas not computed")

    df["rank_in_target"] = df.groupby("target")["rank_metric"].rank(method="min")
    df = df.sort_values(["target", "rank_metric"])

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "docking_report.xlsx"

    summary_cols = [c for c in ["target", "ligand", "rank_in_target",
                                "rank_metric", "delta_vs_reference",
                                "cluster_id", "cluster_size", "cluster_sd",
                                "best_score", "top4_mean", "n_clusters",
                                "n_structures", "engine", "note"]
                    if c in df.columns]

    notes = pd.DataFrame({
        "note": [
            f"Engine: {engine}",
            "rank_metric = mean HADDOCK score of the best 4 members of the "
            "best-scoring cluster; falls back to top4_mean when no cluster "
            "analysis is present.",
            "More negative is better.",
            "delta_vs_reference is computed WITHIN a target only. HADDOCK "
            "scores are not comparable across different targets.",
            "The HADDOCK score is an empirical ranking function, NOT a binding "
            "affinity. It does not give a Kd.",
            "Cluster size is evidence: a large, tight cluster is a more "
            "reproducible binding mode than a small, scattered one.",
            "Scores from HADDOCK 2.5 and HADDOCK3 must not be pooled.",
        ]})

    with pd.ExcelWriter(dest, engine="openpyxl") as writer:
        df[summary_cols].to_excel(writer, sheet_name="Summary", index=False)
        df.to_excel(writer, sheet_name="Full", index=False)
        for tname, grp in df.groupby("target"):
            sheet = f"T_{tname}"[:31]
            grp[summary_cols].to_excel(writer, sheet_name=sheet, index=False)
        notes.to_excel(writer, sheet_name="HowToRead", index=False)

    print()
    for tname, grp in df.groupby("target"):
        print(f"  [{tname}] best 5 by cluster score:")
        for _, r in grp.head(5).iterrows():
            d = r.get("delta_vs_reference")
            dtxt = f"  (Δ vs ref {d:+.2f})" if pd.notna(d) else ""
            print(f"    {r['ligand']:<40} {r['rank_metric']:>10.2f}{dtxt}")
    print()
    print(f"  wrote {dest}")

    write_manifest(out_dir, "score", args.config, [],
                   extra={"engine": engine, "reference": reference,
                          "n_pairs": int(len(df)), "top_n": TOP_N})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
