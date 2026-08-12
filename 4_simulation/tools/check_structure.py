#!/usr/bin/env python3
"""
check_structure.py — pre-flight geometry check for ANY input PDB.

Run this BEFORE spending queue time. It reports the three defect classes that
actually break GROMACS, and nothing else:

  1. D-amino-acid centres  — inverted stereocentres. A minimiser can never fix
     one (inversion is bond-breaking), so these MUST be caught before setup.
  2. Hard steric overlaps  — heavy atoms closer than ~2.0 A. Below ~1.2 A the
     Lennard-Jones term goes to infinity and EM reports Epot = inf.
  3. Ring piercings        — a bond threaded through an aromatic ring. Survives
     minimisation and restrained NVT, then detonates in free dynamics.

Usage:
    python3 check_structure.py file.pdb [file2.pdb ...]

Requires: biopython, numpy   (pip install --user biopython numpy)

INTERPRETING THE OUTPUT
    D-centres 0, piercings 0, severe overlaps 0   -> go straight to setup
    a few contacts 1.8-2.3 A that are N/O pairs   -> normal for model structures;
                                                     EM + the NPT warm-up handle it
    any D-centre or piercing                      -> fix the structure first;
                                                     see docs/STRUCTURE_PREP.md
"""
import sys
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

SEVERE = 1.2     # A — essentially superimposed
CLOSE = 2.0      # A — worth reporting
RINGS = {
    'PHE': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    'TYR': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    'HIS': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
    'TRP': [['CG', 'CD1', 'NE1', 'CE2', 'CD2'],
            ['CD2', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3']],
}


def improper_deg(res):
    """Signed N-CA-C-CB improper. L-amino acids all share one sign."""
    ca = res['CA'].coord
    b1, b2, b3 = res['N'].coord - ca, res['C'].coord - ca, res['CB'].coord - ca
    n1 = np.cross(b1, b2)
    return float(np.degrees(np.arctan2(
        np.dot(np.cross(n1, b3), b2 / np.linalg.norm(b2)), np.dot(n1, b3))))


def check(path):
    s = PDBParser(QUIET=True).get_structure('x', path)[0]

    # --- chirality ---
    rs = [r for r in s.get_residues() if {'N', 'CA', 'C', 'CB'} <= {a.name for a in r}]
    angs = [improper_deg(r) for r in rs]
    med = np.median(angs) if angs else 0.0
    D = [(r.get_parent().id, r.id[1], r.resname)
         for r in rs if np.sign(improper_deg(r)) != np.sign(med)]

    # --- clashes ---
    atoms = [a for a in s.get_atoms() if a.element != 'H']
    ns = NeighborSearch(atoms)
    close, severe = [], []
    for a, b in ns.search_all(CLOSE):
        ra, rb = a.get_parent(), b.get_parent()
        if ra is rb:
            continue
        if ra.get_parent().id == rb.get_parent().id and abs(ra.id[1] - rb.id[1]) == 1:
            continue
        d = float(a - b)
        rec = (round(d, 2), ra.get_parent().id, ra.resname, ra.id[1], a.name,
               rb.get_parent().id, rb.resname, rb.id[1], b.name)
        close.append(rec)
        if d < SEVERE:
            severe.append(rec)
    close.sort()

    # --- ring piercings ---
    pierce = []
    for res in s.get_residues():
        for ring in RINGS.get(res.resname, []):
            try:
                ra = [res[n] for n in ring]
            except KeyError:
                continue
            rc = np.array([x.coord for x in ra])
            ctr = rc.mean(0)
            nrm = np.zeros(3)
            for i in range(len(rc)):
                nrm += np.cross(rc[i] - ctr, rc[(i + 1) % len(rc)] - ctr)
            nrm /= np.linalg.norm(nrm)
            rad = np.max(np.linalg.norm(rc - ctr, axis=1))
            rset = {id(x) for x in ra}
            for a in atoms:
                if id(a) in rset or a.get_parent() is res:
                    continue
                dz = np.dot(a.coord - ctr, nrm)
                if abs(dz) > 1.6:
                    continue
                if np.linalg.norm(a.coord - ctr - dz * nrm) < rad * 0.8:
                    pierce.append((res.resname, res.id[1],
                                   a.get_parent().resname, a.get_parent().id[1], a.name))

    chains = sorted({c.id for c in s})
    print(f"\n=== {path} ===")
    print(f"  chains              : {chains}  ({len(chains)})")
    print(f"  heavy atoms         : {len(atoms)}")
    print(f"  hydrogens present   : {any(a.element == 'H' for a in s.get_atoms())}")
    print(f"  D-centres           : {len(D)}      <- must be 0")
    for x in D[:10]:
        print(f"       {x}")
    print(f"  severe overlaps<{SEVERE}A: {len(severe)}      <- must be 0")
    print(f"  contacts <{CLOSE}A      : {len(close)}      <- a few N/O pairs is normal")
    for c in close[:8]:
        print("       %.2fA %s/%s%d:%s -- %s/%s%d:%s" % c)
    print(f"  ring piercings      : {len(pierce)}      <- must be 0")
    for h in pierce[:5]:
        print(f"       {h}")

    ok = not D and not severe and not pierce
    print(f"  VERDICT             : {'PASS - ready for setup' if ok else 'FIX BEFORE SETUP'}")
    return ok


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    results = [check(p) for p in sys.argv[1:]]
    sys.exit(0 if all(results) else 1)
