#!/usr/bin/env python3
"""
01_geometry_check.py — hard geometry gate.

Checks the three defect classes that actually break downstream docking and MD,
and nothing else. Each is something no energy minimiser can repair, which is why
they must be caught here rather than discovered after a week of compute:

  1. D-amino-acid centres — an inverted stereocentre. Correcting one means
     breaking and re-forming bonds, so minimisation can never fix it. It
     survives energy minimisation and restrained equilibration, then detonates
     the moment dynamics are unrestrained.

  2. Hard steric overlaps — heavy atoms closer than ~1.2 A. The Lennard-Jones
     r^-12 term diverges, and GROMACS reports Epot = inf during minimisation.

  3. Ring piercings — a bond threaded through an aromatic ring. Topologically a
     knot; minimisation relaxes around it rather than undoing it.

Mild contacts (1.8-2.3 A between N/O atoms) are NOT flagged as failures. They
are normal in predicted structures and are resolved by ordinary energy
minimisation plus a gentle equilibration warm-up. Treating them as defects is
what leads people to "fix" structures that were never broken.

Usage:
    python3 01_geometry_check.py config/hadv_c5_hvr7.yaml [--models DIR]

Writes: reports/validation/geometry.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

try:
    from Bio.PDB import PDBParser, NeighborSearch, MMCIFParser
except ImportError:
    sys.exit("Biopython missing.  pip install -r requirements.txt")

from common import (banner, load_config, model_id, write_manifest)

SEVERE_CUTOFF = 1.2   # A — essentially superimposed
CLOSE_CUTOFF = 2.0    # A — reported for information only

RINGS = {
    "PHE": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
    "TYR": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
    "HIS": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "TRP": [["CG", "CD1", "NE1", "CE2", "CD2"],
            ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]],
}


def improper_deg(res) -> float:
    """Signed N-CA-C-CB improper dihedral. All L-residues share one sign."""
    ca = res["CA"].coord
    b1 = res["N"].coord - ca
    b2 = res["C"].coord - ca
    b3 = res["CB"].coord - ca
    n1 = np.cross(b1, b2)
    return float(np.degrees(np.arctan2(
        np.dot(np.cross(n1, b3), b2 / np.linalg.norm(b2)), np.dot(n1, b3))))


def load_structure(path: Path):
    parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
    return parser.get_structure("m", str(path))[0]


def check_chirality(model):
    residues = [r for r in model.get_residues()
                if {"N", "CA", "C", "CB"} <= {a.name for a in r}]
    if not residues:
        return [], 0.0
    angles = [improper_deg(r) for r in residues]
    median = float(np.median(angles))
    flagged = [(r.get_parent().id, r.id[1], r.resname, round(improper_deg(r), 1))
               for r in residues if np.sign(improper_deg(r)) != np.sign(median)]
    return flagged, median


def check_clashes(model):
    atoms = [a for a in model.get_atoms() if a.element != "H"]
    ns = NeighborSearch(atoms)
    close, severe = [], []
    for a, b in ns.search_all(CLOSE_CUTOFF):
        ra, rb = a.get_parent(), b.get_parent()
        if ra is rb:
            continue
        # sequential residues in the same chain are bonded, not clashing
        if ra.get_parent().id == rb.get_parent().id and abs(ra.id[1] - rb.id[1]) == 1:
            continue
        d = float(a - b)
        rec = (round(d, 2),
               f"{ra.get_parent().id}/{ra.resname}{ra.id[1]}:{a.name}",
               f"{rb.get_parent().id}/{rb.resname}{rb.id[1]}:{b.name}")
        close.append(rec)
        if d < SEVERE_CUTOFF:
            severe.append(rec)
    close.sort()
    return close, severe, len(atoms)


def check_piercings(model):
    atoms = [a for a in model.get_atoms() if a.element != "H"]
    hits = []
    for res in model.get_residues():
        for ring in RINGS.get(res.resname, []):
            try:
                ring_atoms = [res[n] for n in ring]
            except KeyError:
                continue
            coords = np.array([x.coord for x in ring_atoms])
            centre = coords.mean(0)
            normal = np.zeros(3)
            for i in range(len(coords)):
                normal += np.cross(coords[i] - centre,
                                   coords[(i + 1) % len(coords)] - centre)
            norm = np.linalg.norm(normal)
            if norm == 0:
                continue
            normal /= norm
            radius = float(np.max(np.linalg.norm(coords - centre, axis=1)))
            ring_ids = {id(x) for x in ring_atoms}
            for a in atoms:
                if id(a) in ring_ids or a.get_parent() is res:
                    continue
                dz = float(np.dot(a.coord - centre, normal))
                if abs(dz) > 1.6:
                    continue
                if np.linalg.norm(a.coord - centre - dz * normal) < radius * 0.8:
                    hits.append(
                        f"{res.resname}{res.id[1]} pierced by "
                        f"{a.get_parent().resname}{a.get_parent().id[1]}:{a.name}")
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--models", default=None,
                    help="directory of model PDB/CIF files "
                         "(default: paths.models from config)")
    ap.add_argument("--out", default="reports/validation")
    args = ap.parse_args()

    cfg = load_config(args.config)
    thr = cfg["validation"]["thresholds"]

    model_dir = Path(args.models or (cfg.get("paths") or {}).get("models", "data/models"))
    if not model_dir.is_dir():
        sys.exit(f"[FATAL] model directory not found: {model_dir}\n"
                 f"        pass --models DIR or set paths.models in the config")

    files = sorted([p for p in model_dir.iterdir()
                    if p.suffix in (".pdb", ".cif")])
    if not files:
        sys.exit(f"[FATAL] no .pdb/.cif files in {model_dir}")

    banner(f"Geometry check — {len(files)} model(s) from {model_dir}")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []

    for path in files:
        mid = model_id(path)
        try:
            model = load_structure(path)
        except Exception as exc:
            print(f"  [ERROR] {mid}: {exc}")
            rows.append({"model": mid, "geometry_pass": False,
                         "error": str(exc)[:120]})
            continue

        chains = sorted({c.id for c in model})
        d_flagged, _ = check_chirality(model)
        close, severe, n_atoms = check_clashes(model)
        piercings = check_piercings(model)

        checks = {
            "d_centres": len(d_flagged) <= thr["max_d_centres"],
            "severe_clashes": len(severe) <= thr["max_severe_clashes"],
            "piercings": len(piercings) <= thr["max_ring_piercings"],
            "chains": len(chains) == thr["expected_chains"],
        }
        passed = all(checks.values())
        failed = [k for k, v in checks.items() if not v]

        rows.append({
            "model": mid,
            "n_chains": len(chains),
            "chains": ",".join(chains),
            "heavy_atoms": n_atoms,
            "d_centres": len(d_flagged),
            "severe_clashes": len(severe),
            "close_contacts": len(close),
            "ring_piercings": len(piercings),
            "geometry_pass": passed,
            "geometry_fail_reason": ";".join(failed),
            "worst_contact_A": close[0][0] if close else "",
            "d_centre_detail": ";".join(f"{c}/{r}{n}" for c, n, r, _ in d_flagged[:6]),
            "piercing_detail": ";".join(piercings[:3]),
        })

        flag = "PASS" if passed else "FAIL"
        print(f"  [{flag}] {mid:<40} chains={len(chains)} "
              f"D={len(d_flagged)} severe={len(severe)} pierce={len(piercings)}")
        if not passed:
            print(f"         reason: {', '.join(failed)}")

    dest = out_dir / "geometry.csv"
    with open(dest, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    n_pass = sum(1 for r in rows if r.get("geometry_pass"))
    print()
    print(f"  {n_pass}/{len(rows)} passed geometry")
    print(f"  wrote {dest}")

    write_manifest(out_dir, "geometry", args.config, files,
                   extra={"n_pass": n_pass, "n_total": len(rows)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
