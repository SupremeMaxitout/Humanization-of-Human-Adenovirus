#!/usr/bin/env python3
"""
02_homology_search.py — find human sequence windows homologous to a viral region.

Slides a window of `homology.window_size` across every human protein, scores it
against the wild-type region, and keeps windows that clear both an alignment
threshold and a physicochemical (charge) filter.

Everything tunable lives in the config:
    window_size, substitution_matrix, alignment_mode, gap scores,
    score_threshold, ph, charge_deviation, match_polarity, dedupe_on

CHANGES vs the original script01.py
  * compare_with_QUERY was being called twice per protein (once inside the
    search, once in main). Now called once.
  * A full pairwise alignment was computed and printed for every protein purely
    for display, then discarded. Removed — that was the dominant cost across
    ~110k proteins.
  * alignment_mode now defaults to 'local'. Globally aligning a 15-residue
    window against a 41-residue query penalises the length difference as
    end gaps, so scores mostly reflected length, not similarity. Local
    alignment compares like with like. Set 'global' in config to reproduce
    the original behaviour.
  * Reports normalised score (score/window_size) so thresholds transfer
    across different window sizes.
  * Progress is periodic rather than a wall of text.

Usage:
    python3 02_homology_search.py config/hadv_c5_hvr7.yaml
    python3 02_homology_search.py config/hadv_c5_hvr7.yaml --limit 500
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, write_manifest, banner, ensure_dirs

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def make_aligner(hcfg: dict) -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load(hcfg["substitution_matrix"])
    al.mode = hcfg.get("alignment_mode", "local")
    al.open_gap_score = hcfg.get("open_gap_score", -10)
    al.extend_gap_score = hcfg.get("extend_gap_score", -0.5)
    return al


def charge_profile(seq: str, ph: float) -> dict:
    pa = ProteinAnalysis(seq)
    return {"pI": round(pa.isoelectric_point(), 2),
            "net_charge": round(pa.charge_at_pH(ph), 4)}


def best_window(seq: str, query: str, aligner, window: int):
    """Highest-scoring window of length `window` in `seq` against `query`."""
    best_score = -float("inf")
    best_sub = ""
    n = len(seq) - window + 1
    for i in range(n):
        sub = seq[i:i + window]
        if not set(sub) <= VALID_AA:      # skip X, U, * and other odd codes
            continue
        try:
            score = aligner.score(sub, query)
        except Exception:
            continue
        if score > best_score:
            best_score = score
            best_sub = sub
    return best_sub, best_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    h = cfg["homology"]
    query = cfg["region"]["query"].strip().upper()
    ensure_dirs(cfg, "reports")

    proteome_path = Path(h["proteome_json"])
    if not proteome_path.exists():
        sys.exit(f"[ERROR] {proteome_path} not found. Run 01_parse_proteome.py first.")

    aligner = make_aligner(h)
    q = charge_profile(query, h["ph"])
    q_sign = 1 if q["net_charge"] > 0 else (-1 if q["net_charge"] < 0 else 0)

    banner(f"Homology search — {cfg['region']['name']} of {cfg['target']['virus']}")
    print(f"  query           : {query}  ({len(query)} aa)")
    print(f"  net charge @pH{h['ph']} : {q['net_charge']}   pI {q['pI']}")
    print(f"  window          : {h['window_size']} aa")
    print(f"  matrix / mode   : {h['substitution_matrix']} / {aligner.mode}")
    print(f"  score threshold : {h['score_threshold']}")
    print(f"  charge tol      : +/- {h['charge_deviation']}")
    print(f"  polarity filter : {h['match_polarity']}\n")

    with open(proteome_path, encoding="utf-8") as fh:
        proteome = json.load(fh)

    items = list(proteome.items())
    if args.limit:
        items = items[:args.limit]

    leads = []
    n_scanned = n_score_pass = n_polarity_fail = n_charge_fail = 0

    for i, (ensp, entry) in enumerate(items, 1):
        seq = (entry.get("sequence") or "").strip().upper()
        if len(seq) < h["window_size"]:
            continue
        n_scanned += 1

        if i % 5000 == 0:
            print(f"  ... {i}/{len(items)} scanned, {len(leads)} hits")

        sub, score = best_window(seq, query, aligner, h["window_size"])
        if not sub or score < h["score_threshold"]:
            continue
        n_score_pass += 1

        c = charge_profile(sub, h["ph"])
        if h.get("match_polarity", True) and q_sign != 0:
            c_sign = 1 if c["net_charge"] > 0 else (-1 if c["net_charge"] < 0 else 0)
            if c_sign != q_sign:
                n_polarity_fail += 1
                continue

        diff = abs(q["net_charge"] - c["net_charge"])
        if diff > h["charge_deviation"]:
            n_charge_fail += 1
            continue

        leads.append({
            "protein_ID": ensp,
            "description": entry.get("description", ""),
            "candidate_sequence": sub,
            "score": round(float(score), 2),
            "score_per_residue": round(float(score) / h["window_size"], 3),
            "net_charge": c["net_charge"],
            "pI": c["pI"],
            "charge_difference": round(diff, 4),
        })

    print(f"\n  scanned            : {n_scanned}")
    print(f"  passed score       : {n_score_pass}")
    print(f"  failed polarity    : {n_polarity_fail}")
    print(f"  failed charge tol  : {n_charge_fail}")
    print(f"  HITS               : {len(leads)}")

    if not leads:
        print("\n  No candidates. Loosen score_threshold or charge_deviation.")
        return

    df = pd.DataFrame(leads).sort_values("score", ascending=False)
    before = len(df)
    df = df.drop_duplicates(subset=h.get("dedupe_on", "description"), keep="first")
    print(f"  after dedupe on '{h.get('dedupe_on','description')}': "
          f"{len(df)} (removed {before - len(df)})")

    df = df.reset_index(drop=True)
    out = Path(cfg["paths"]["reports"]) / f"{cfg['project']['name']}_leads.xlsx"
    df.to_excel(out, index=False)
    print(f"\n  written: {out}")

    write_manifest(cfg, "02_homology_search",
                   inputs={"proteome_json": str(proteome_path),
                           "n_proteins_scanned": n_scanned},
                   outputs={"xlsx": str(out), "n_leads": len(df)},
                   extra={"query": query, "query_charge": q,
                          "window_size": h["window_size"],
                          "matrix": h["substitution_matrix"],
                          "mode": aligner.mode,
                          "score_threshold": h["score_threshold"]})


if __name__ == "__main__":
    main()
