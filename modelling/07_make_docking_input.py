#!/usr/bin/env python3
"""
07_make_docking_input.py — merge the homotrimer into ONE chain for docking.

HADDOCK treats each chain as a separate body. A homotrimer receptor should be
presented as a single rigid body, so chains B and C are renumbered by
1x and 2x protomer_length and folded into chain A.

  chain A: residues     1 ..  N
  chain B: residues N+1 .. 2N   (was 1..N)
  chain C: residues 2N+1 .. 3N  (was 1..N)

*** READ THIS ***
These single-chain files are for DOCKING ONLY.

Never feed them to GROMACS. pdb2gmx will create real peptide bonds across the
protomer junctions, artificially stiffening the exact inter-protomer interface
an HVR flexibility study measures. That mistake produced three unusable MD
trajectories on this project. MD uses the multi-chain files from step 06.

Shift size comes from target.protomer_length in the config rather than being
hardcoded, so it follows the target.

Usage:
    python3 07_make_docking_input.py config/hadv_c5_hvr7.yaml
"""
import argparse
import sys
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, write_manifest, banner, ensure_dirs

ap = argparse.ArgumentParser()
ap.add_argument("config", nargs="?", default="config/hadv_c5_hvr7.yaml")
args = ap.parse_args()

cfg = load_config(args.config)
ensure_dirs(cfg, "pdb_singlechain")

src_dir = Path(cfg["paths"]["pdb_multichain"])
dst_dir = Path(cfg["paths"]["pdb_singlechain"])
plen = int(cfg["target"]["protomer_length"])
nch = int(cfg["target"]["n_chains"])

# chain letter -> residue offset, derived from protomer length
letters = [chr(ord("A") + i) for i in range(nch)]
SHIFTS = {letters[i]: plen * i for i in range(1, nch)}

banner("Building single-chain docking inputs")
print(f"  protomer length : {plen}")
print(f"  shifts          : {SHIFTS}\n")

parser = PDBParser(QUIET=True)
io = PDBIO()
ok, failed = 0, []

for pdb in sorted(src_dir.glob("*.pdb")):
    try:
        st = parser.get_structure(pdb.stem, str(pdb))
        model = next(st.get_models())
        present = [c.id for c in model]

        if "A" not in present:
            raise ValueError("no chain A to merge into")

        moved = []
        for letter, shift in SHIFTS.items():
            if letter not in present:
                continue
            ch = model[letter]
            for res in list(ch):
                het, num, icode = res.id
                ch.detach_child(res.id)
                res.id = (het, num + shift, icode)
                moved.append(res)

        chain_a = model["A"]
        for res in moved:
            chain_a.add(res)
        for letter in SHIFTS:
            if letter in [c.id for c in model]:
                model.detach_child(letter)

        chain_a.child_list.sort(key=lambda r: r.id[1])

        out = dst_dir / f"{pdb.stem}_docking.pdb"
        io.set_structure(st)
        io.save(str(out))

        nres = len(chain_a.child_list)
        print(f"  [ok] {pdb.name} -> {out.name}  ({nres} residues in chain A)")
        if nres != plen * nch:
            print(f"       [WARN] expected {plen * nch} residues, got {nres}")
        ok += 1
    except Exception as e:
        print(f"  [FAIL] {pdb.name}: {e}")
        failed.append(pdb.name)

print(f"\n  written : {ok}   failed: {len(failed)}")
print(f"  output  : {dst_dir}")
print("\n  REMINDER: single-chain = docking only. MD uses the multi-chain set.")

write_manifest(cfg, "07_make_docking_input",
               inputs={"pdb_multichain": str(src_dir)},
               outputs={"pdb_singlechain": str(dst_dir), "n_written": ok},
               extra={"shifts": SHIFTS, "failed": failed})
