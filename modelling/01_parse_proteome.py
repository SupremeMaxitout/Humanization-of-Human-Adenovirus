#!/usr/bin/env python3
"""
01_parse_proteome.py — Ensembl peptide FASTA -> JSON keyed by ENSP.

Fixes two bugs from the original:
  * json.dump(..., indecnt=2)  ->  indent=2   (crashed on write)
  * the `limit` early-break skipped the final record

Get the proteome first (~15 MB gz, ~110k sequences):
    mkdir -p data && cd data
    wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
    gunzip Homo_sapiens.GRCh38.pep.all.fa.gz

Usage:
    python3 01_parse_proteome.py config/hadv_c5_hvr7.yaml [--limit N]
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, write_manifest, banner, file_sha256

ap = argparse.ArgumentParser()
ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
ap.add_argument("--limit", type=int, default=None,
                help="parse only the first N proteins (for a quick test)")
args = ap.parse_args()

cfg = load_config(args.config)
fasta = Path(cfg["homology"]["proteome_fasta"])
out_json = Path(cfg["homology"]["proteome_json"])

if not fasta.exists():
    sys.exit(f"[ERROR] proteome FASTA not found: {fasta}\n"
             f"        Download it (see this script's docstring).")

banner("Parsing proteome FASTA")

proteins: dict[str, dict] = {}
ensp = desc = None
buf: list[str] = []


def flush():
    if ensp is not None:
        proteins[ensp] = {"description": desc, "sequence": "".join(buf)}


with open(fasta, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            if args.limit and len(proteins) >= args.limit:
                ensp = None
                break
            buf = []
            header = line[1:]
            ensp = header.split()[0].split(".")[0]
            desc = ""
            if "description:" in header:
                desc = header.split("description:", 1)[1].strip()
        else:
            buf.append(line)
    flush()          # <-- the original dropped the last record

out_json.parent.mkdir(parents=True, exist_ok=True)
with open(out_json, "w", encoding="utf-8") as fh:
    json.dump(proteins, fh, indent=2)      # <-- was 'indecnt'

lengths = [len(v["sequence"]) for v in proteins.values()]
print(f"  proteins parsed : {len(proteins)}")
print(f"  median length   : {sorted(lengths)[len(lengths)//2] if lengths else 0} aa")
print(f"  written         : {out_json}")

write_manifest(cfg, "01_parse_proteome",
               inputs={"fasta": str(fasta), "sha256": file_sha256(fasta)},
               outputs={"json": str(out_json), "n_proteins": len(proteins)})
