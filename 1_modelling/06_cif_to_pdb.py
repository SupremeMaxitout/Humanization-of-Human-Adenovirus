#!/usr/bin/env python3
"""
06_cif_to_pdb.py — mmCIF -> PDB, preserving chains and pLDDT.

Fixes from the original:
  * output was named "<file>.cif.pdb" because the .cif extension was not
    stripped. Names containing '.cif' mid-string make PDBFixer and other
    format-sniffing tools try to parse a PDB as mmCIF, which fails with
    "IndexError: list index out of range". Now the suffix is replaced properly.
  * uses gemmi.read_structure() rather than read_file + make_structure_from_block
  * setup_entities() + assign_label_seq_id() so chain records survive
  * ColabFold output is already PDB, so it is copied through unchanged

pLDDT lives in the B-factor column and is carried over — Tool 2 reads it there.

Usage:
    python3 06_cif_to_pdb.py config/hadv_c5_hvr7.yaml
"""
import argparse
import shutil
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, write_manifest, banner, ensure_dirs

ap = argparse.ArgumentParser()
ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
args = ap.parse_args()

cfg = load_config(args.config)
ensure_dirs(cfg, "pdb_multichain")

src_dir = Path(cfg["paths"]["best_models"])
dst_dir = Path(cfg["paths"]["pdb_multichain"])
expect_chains = int(cfg["target"]["n_chains"])

banner("Converting best models to PDB (multi-chain, for MD)")

ok, failed = 0, []
for f in sorted(src_dir.iterdir()):
    if f.suffix.lower() not in (".cif", ".pdb"):
        continue
    out = dst_dir / (f.stem + ".pdb")     # <- .cif properly stripped
    try:
        if f.suffix.lower() == ".pdb":
            shutil.copy2(f, out)
            st = gemmi.read_structure(str(out))
        else:
            st = gemmi.read_structure(str(f))
            st.setup_entities()
            st.assign_label_seq_id()
            st.write_pdb(str(out))

        chains = [c.name for c in st[0]]
        flag = "" if len(chains) == expect_chains else \
               f"  [WARN] expected {expect_chains} chains, got {len(chains)}"
        print(f"  [ok] {f.name} -> {out.name}  chains={','.join(chains)}{flag}")
        ok += 1
    except Exception as e:
        print(f"  [FAIL] {f.name}: {e}")
        failed.append(f.name)

print(f"\n  converted : {ok}\n  failed    : {len(failed)}")
print(f"  output    : {dst_dir}")
print("\n  These MULTI-CHAIN files are for MD. Docking uses the single-chain")
print("  files from 07_make_docking_input.py — do not mix them up.")

write_manifest(cfg, "06_cif_to_pdb",
               inputs={"best_models": str(src_dir)},
               outputs={"pdb_multichain": str(dst_dir), "n_converted": ok},
               extra={"failed": failed})
