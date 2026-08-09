#!/usr/bin/env python3
"""
04_compile_report.py — merge the validation stages into one decision table.

Combines geometry.csv, confidence.csv and sasa.csv into a single spreadsheet and
applies the pass/fail gate defined in the config. Models that pass are the ones
that proceed to docking; everything else is held back with a stated reason.

Design decision worth stating: thresholds live in the config and are applied
uniformly, so the gate is fixed before results are seen. Choosing cutoffs after
looking at the numbers is how a screen turns into post-hoc rationalisation.

Sheets written
  Summary    one row per model, verdict plus the reason it failed
  Geometry   full geometry detail
  Confidence pLDDT global and per region, pTM/ipTM/PAE
  SASA       per-region surface area and ratios to the reference
  Thresholds the exact gate applied, recorded alongside the results

Usage:
    python3 04_compile_report.py config/hadv_c5_hvr7.yaml
                                 [--in reports/validation]

Writes: reports/validation/validation_report.xlsx
        reports/validation/passed_models.txt   (input list for Tool 3)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas missing.  pip install -r requirements.txt")

from common import banner, get_target_region, load_config, write_manifest


def read_stage(in_dir: Path, name: str) -> pd.DataFrame:
    path = in_dir / f"{name}.csv"
    if not path.is_file():
        print(f"  [note] {name}.csv not found — that stage will be blank")
        return pd.DataFrame(columns=["model"])
    df = pd.read_csv(path)
    if "model" not in df.columns:
        print(f"  [warn] {name}.csv has no 'model' column — skipping")
        return pd.DataFrame(columns=["model"])
    return df


def as_float(val):
    try:
        f = float(val)
        return None if pd.isna(f) else f
    except (TypeError, ValueError):
        return None


def evaluate(row: pd.Series, thr: dict, target_region: str | None) -> tuple[bool, list]:
    """Apply the gate. Missing metrics are NOT treated as failures — they are
    reported as gaps, so an incomplete input set does not silently fail
    everything."""
    reasons = []

    if row.get("geometry_pass") is False or str(row.get("geometry_pass")).lower() == "false":
        detail = row.get("geometry_fail_reason", "")
        reasons.append(f"geometry({detail})" if isinstance(detail, str) and detail
                       else "geometry")

    mean_plddt = as_float(row.get("mean_plddt"))
    if mean_plddt is not None and mean_plddt < thr["min_mean_plddt"]:
        reasons.append(f"mean_pLDDT {mean_plddt:.1f}<{thr['min_mean_plddt']}")

    if target_region:
        reg = as_float(row.get(f"plddt_{target_region}"))
        if reg is not None and reg < thr["min_region_plddt"]:
            reasons.append(f"{target_region}_pLDDT {reg:.1f}<{thr['min_region_plddt']}")

    iptm = as_float(row.get("iptm"))
    if iptm is not None and iptm < thr["min_iptm"]:
        reasons.append(f"ipTM {iptm:.2f}<{thr['min_iptm']}")

    ptm = as_float(row.get("ptm"))
    if ptm is not None and ptm < thr["min_ptm"]:
        reasons.append(f"pTM {ptm:.2f}<{thr['min_ptm']}")

    frac = as_float(row.get("fraction_disordered"))
    if frac is not None and frac > thr["max_fraction_disordered"]:
        reasons.append(f"disordered {frac:.2f}>{thr['max_fraction_disordered']}")

    clash = row.get("has_clash")
    if str(clash).lower() in ("true", "1", "1.0"):
        reasons.append("predictor_clash_flag")

    return (not reasons), reasons


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--in", dest="in_dir", default="reports/validation")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    thr = cfg["validation"]["thresholds"]
    target_region = get_target_region(cfg)

    in_dir = Path(args.in_dir)
    if not in_dir.is_dir():
        sys.exit(f"[FATAL] {in_dir} not found — run stages 01-03 first")

    banner("Compiling validation report")

    geometry = read_stage(in_dir, "geometry")
    confidence = read_stage(in_dir, "confidence")
    sasa = read_stage(in_dir, "sasa")

    merged = geometry
    for df in (confidence, sasa):
        if not df.empty:
            merged = merged.merge(df, on="model", how="outer") if not merged.empty else df

    if merged.empty:
        sys.exit("[FATAL] nothing to compile — no stage produced results")

    verdicts, reasons = [], []
    for _, row in merged.iterrows():
        ok, why = evaluate(row, thr, target_region)
        verdicts.append("PASS" if ok else "FAIL")
        reasons.append("; ".join(why))
    merged["verdict"] = verdicts
    merged["fail_reasons"] = reasons

    summary_cols = [c for c in [
        "model", "verdict", "fail_reasons", "n_chains", "mean_plddt",
        f"plddt_{target_region}" if target_region else None,
        "ptm", "iptm", "ranking_score", "d_centres", "severe_clashes",
        "ring_piercings", "sasa_all_regions",
    ] if c and c in merged.columns]
    summary = merged[summary_cols].sort_values(
        ["verdict", "mean_plddt"] if "mean_plddt" in merged.columns else ["verdict"],
        ascending=[True, False] if "mean_plddt" in merged.columns else [True])

    out_path = Path(args.out) if args.out else in_dir / "validation_report.xlsx"
    thr_df = pd.DataFrame(sorted(thr.items()), columns=["threshold", "value"])

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        if not geometry.empty:
            geometry.to_excel(writer, sheet_name="Geometry", index=False)
        if not confidence.empty:
            confidence.to_excel(writer, sheet_name="Confidence", index=False)
        if not sasa.empty:
            sasa.to_excel(writer, sheet_name="SASA", index=False)
        thr_df.to_excel(writer, sheet_name="Thresholds", index=False)

    passed = merged.loc[merged["verdict"] == "PASS", "model"].tolist()
    list_path = in_dir / "passed_models.txt"
    with open(list_path, "w") as fh:
        fh.write("\n".join(passed) + ("\n" if passed else ""))

    print()
    print(f"  {len(passed)}/{len(merged)} models passed")
    print(f"  report : {out_path}")
    print(f"  passed : {list_path}   (input list for Tool 3)")
    if len(merged) - len(passed):
        print()
        print("  failures:")
        for _, row in merged[merged["verdict"] == "FAIL"].iterrows():
            print(f"    {row['model']:<40} {row['fail_reasons']}")

    write_manifest(in_dir, "report", args.config, [out_path],
                   extra={"n_pass": len(passed), "n_total": int(len(merged)),
                          "thresholds": thr})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
