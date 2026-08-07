#!/usr/bin/env python3
"""
03_make_prediction_inputs.py — leads spreadsheet -> AF3 JSON and/or ColabFold FASTA.

Writes BOTH formats so the run is not tied to one predictor:

  * AlphaFold3 server : af3_batch_upload.json  (manual upload — alphafoldserver.com
                        has no public submission API, and its terms prohibit
                        scripted access)
  * ColabFold         : candidates.fasta       (colabfold_batch, fully scriptable)

Homotrimers: AF3 uses "count": 3; ColabFold repeats the sequence 3x separated
by ':'.

Usage:
    python3 03_make_prediction_inputs.py config/hadv_c5_hvr7.yaml
    python3 03_make_prediction_inputs.py config/hadv_c5_hvr7.yaml --leads my.xlsx
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import (load_config, load_template, build_sequence, candidate_id,
                    write_manifest, banner, ensure_dirs)

ap = argparse.ArgumentParser()
ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
ap.add_argument("--leads", default=None, help="override leads .xlsx path")
ap.add_argument("--sheet", default=0, help="sheet name or index")
args = ap.parse_args()

cfg = load_config(args.config)
ensure_dirs(cfg, "reports")

leads_path = Path(args.leads) if args.leads else \
    Path(cfg["paths"]["reports"]) / f"{cfg['project']['name']}_leads.xlsx"
if not leads_path.exists():
    sys.exit(f"[ERROR] leads file not found: {leads_path}\n"
             f"        Run 02_homology_search.py first, or pass --leads")

sheet = args.sheet
try:
    sheet = int(sheet)
except (TypeError, ValueError):
    pass

df = pd.read_excel(leads_path, sheet_name=sheet)
template = load_template(cfg)
copies = int(cfg["prediction"]["copies"])

banner(f"Building prediction inputs — {len(df)} candidates")

af3_jobs = []
fasta_lines = []
index_rows = []

for i, row in df.iterrows():
    graft = str(row["candidate_sequence"]).strip().upper()
    seq = build_sequence(template, graft)
    name = candidate_id(i, row["protein_ID"])

    af3_jobs.append({
        "name": name,
        "modelSeeds": [],
        "sequences": [{"proteinChain": {"sequence": seq, "count": copies}}],
    })

    fasta_lines.append(f">{name}")
    fasta_lines.append(":".join([seq] * copies))

    index_rows.append({
        "candidate": name,
        "protein_ID": row["protein_ID"],
        "description": row.get("description", ""),
        "candidate_sequence": graft,
        "full_length": len(seq),
        "score": row.get("score"),
    })

reports = Path(cfg["paths"]["reports"])
af3_out = reports / "af3_batch_upload.json"
fasta_out = reports / "candidates.fasta"
index_out = reports / "candidate_index.csv"

af3_out.write_text(json.dumps(af3_jobs, indent=4))
fasta_out.write_text("\n".join(fasta_lines) + "\n")
pd.DataFrame(index_rows).to_csv(index_out, index=False)

lengths = {len(j["sequences"][0]["proteinChain"]["sequence"]) for j in af3_jobs}
print(f"  candidates        : {len(af3_jobs)}")
print(f"  copies per job    : {copies}")
print(f"  full-length range : {min(lengths)}-{max(lengths)} aa")
if len(lengths) > 1:
    print("  [NOTE] candidates differ in length — expected if grafts vary,")
    print("         but it means residue numbering shifts between models.")
print(f"\n  AF3 JSON      : {af3_out}")
print(f"  ColabFold FASTA: {fasta_out}")
print(f"  index (keeps candidate -> source protein mapping): {index_out}")
print("\n  AF3 : upload the JSON at https://alphafoldserver.com (manual)")
print("  CF  : qsub the colabfold job with candidates.fasta as input")

write_manifest(cfg, "03_make_prediction_inputs",
               inputs={"leads": str(leads_path), "n_candidates": len(df)},
               outputs={"af3_json": str(af3_out), "fasta": str(fasta_out),
                        "index": str(index_out)},
               extra={"copies": copies, "template": cfg["target"]["template"]})
