#!/usr/bin/env python3
"""
04_unpack_predictions.py — unpack predictor output into a uniform layout.

Handles BOTH engines:
  * AlphaFold3 server : one .zip per job -> extracted/<job>/
  * ColabFold         : a flat results directory -> normalised into the same shape

Replaces the original hardcoded `for i in range(74)` loop, which silently
skipped anything outside that range and assumed a fixed job count.

Usage:
    python3 04_unpack_predictions.py config/hadv_c5_hvr7.yaml
    python3 04_unpack_predictions.py config/hadv_c5_hvr7.yaml --engine colabfold
"""
import argparse
import shutil
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, write_manifest, banner, ensure_dirs

ap = argparse.ArgumentParser()
ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
ap.add_argument("--engine", default=None, choices=["alphafold3", "colabfold"])
args = ap.parse_args()

cfg = load_config(args.config)
engine = args.engine or cfg["prediction"]["engine"]
ensure_dirs(cfg, "raw_predictions", "extracted")

raw = Path(cfg["paths"]["raw_predictions"])
out = Path(cfg["paths"]["extracted"])

banner(f"Unpacking predictions — engine: {engine}")
print(f"  from: {raw}\n  to  : {out}\n")

ok, failed = 0, []

if engine == "alphafold3":
    zips = sorted(raw.glob("*.zip"))
    if not zips:
        sys.exit(f"[ERROR] no .zip files in {raw}\n"
                 f"        Download AF3 results there first.")
    for z in zips:
        try:
            dest = out / z.stem
            dest.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(z) as zf:
                zf.extractall(dest)
            print(f"  [ok] {z.name}")
            ok += 1
        except zipfile.BadZipFile:
            print(f"  [BAD ZIP] {z.name}")
            failed.append(z.name)
        except Exception as e:
            print(f"  [FAIL] {z.name}: {e}")
            failed.append(z.name)

else:  # colabfold
    # colabfold_batch writes a flat dir: <name>_unrelaxed_rank_001_*.pdb,
    # <name>_scores_rank_001_*.json, <name>_predicted_aligned_error_v1.json ...
    pdbs = sorted(raw.glob("*.pdb"))
    if not pdbs:
        sys.exit(f"[ERROR] no .pdb files in {raw}\n"
                 f"        Point --raw at the colabfold_batch output directory.")
    names = sorted({p.name.split("_unrelaxed")[0].split("_relaxed")[0]
                    for p in pdbs})
    for name in names:
        dest = out / name
        dest.mkdir(parents=True, exist_ok=True)
        moved = 0
        for f in raw.glob(f"{name}*"):
            if f.is_file():
                shutil.copy2(f, dest / f.name)
                moved += 1
        if moved:
            print(f"  [ok] {name}  ({moved} files)")
            ok += 1
        else:
            failed.append(name)

print(f"\n  unpacked : {ok}")
print(f"  failed   : {len(failed)}")
for f in failed:
    print(f"    - {f}")

write_manifest(cfg, "04_unpack_predictions",
               inputs={"raw_dir": str(raw), "engine": engine},
               outputs={"extracted_dir": str(out), "n_unpacked": ok},
               extra={"failed": failed})
