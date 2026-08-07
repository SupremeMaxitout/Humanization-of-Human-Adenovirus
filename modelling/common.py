"""
common.py — shared config loading, paths, naming, and provenance.

Every stage script imports from here so that region boundaries, protomer
length, and naming conventions have exactly ONE definition.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
def load_config(path: str | Path) -> dict:
    """Load a run config and resolve all paths relative to the repo root."""
    path = Path(path)
    if not path.is_absolute():
        path = REPO_ROOT / path
    with open(path) as fh:
        cfg = yaml.safe_load(fh)

    cfg["_config_path"] = str(path)

    for key, val in cfg.get("paths", {}).items():
        cfg["paths"][key] = str(REPO_ROOT / val)
    for key in ("proteome_fasta", "proteome_json"):
        if key in cfg.get("homology", {}):
            cfg["homology"][key] = str(REPO_ROOT / cfg["homology"][key])
    if "template" in cfg.get("target", {}):
        cfg["target"]["template"] = str(REPO_ROOT / cfg["target"]["template"])

    _validate(cfg)
    return cfg


def _validate(cfg: dict) -> None:
    """Fail loudly on the config mistakes that silently corrupt results."""
    problems = []

    region = cfg.get("region", {})
    regions = cfg.get("regions", {})
    name = region.get("name")

    if name and name in regions:
        if list(regions[name]) != list(region.get("range", [])):
            problems.append(
                f"region.range {region.get('range')} disagrees with "
                f"regions.{name} {regions[name]}. These must match or "
                f"residue numbering will shift downstream."
            )

    for rname, rng in regions.items():
        if len(rng) != 2 or rng[0] > rng[1]:
            problems.append(f"regions.{rname} = {rng} is not a valid [start, end]")

    plen = cfg.get("target", {}).get("protomer_length")
    if plen:
        for rname, rng in regions.items():
            if rng[1] > plen:
                problems.append(
                    f"regions.{rname} ends at {rng[1]} but protomer_length "
                    f"is {plen}"
                )

    if problems:
        print("[CONFIG ERROR]", file=sys.stderr)
        for p in problems:
            print("  - " + p, file=sys.stderr)
        sys.exit(2)


def ensure_dirs(cfg: dict, *keys: str) -> None:
    for k in keys:
        Path(cfg["paths"][k]).mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Template
# --------------------------------------------------------------------------
def load_template(cfg: dict) -> str:
    """Read the scaffold template and verify it has exactly one {} slot."""
    path = Path(cfg["target"]["template"])
    if not path.exists():
        sys.exit(f"[ERROR] template not found: {path}\n"
                 f"        See config/templates/README.md")
    text = path.read_text().strip().replace("\n", "").replace(" ", "")
    n = text.count("{}")
    if n != 1:
        sys.exit(f"[ERROR] template must contain exactly one '{{}}' "
                 f"placeholder, found {n}")
    return text


def build_sequence(template: str, graft: str) -> str:
    return template.format(graft)


# --------------------------------------------------------------------------
# Naming — one convention everywhere
# --------------------------------------------------------------------------
def candidate_id(index: int, protein_id: str) -> str:
    """
    Canonical name for a candidate: 'fold_007_ENSP00000371221'.

    Short and free of dots. Filenames containing '.cif' mid-string make
    format-sniffing tools (PDBFixer, some parsers) misread a PDB as mmCIF.
    """
    clean = str(protein_id).split(".")[0].replace(" ", "_")
    return f"fold_{int(index):03d}_{clean}"


def short_name(candidate: str) -> str:
    """'fold_007_ENSP00000371221' -> 'fold_007' for compact plot labels."""
    parts = candidate.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else candidate


# --------------------------------------------------------------------------
# Provenance — what makes "reproducible" real rather than aspirational
# --------------------------------------------------------------------------
def file_sha256(path: str | Path, limit_mb: int = 500) -> str | None:
    path = Path(path)
    if not path.exists():
        return None
    if path.stat().st_size > limit_mb * 1024 * 1024:
        return "skipped-too-large"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return None


def write_manifest(cfg: dict, stage: str, inputs: dict, outputs: dict,
                   extra: dict | None = None) -> Path:
    """Record what ran, on what, with which versions."""
    reports = Path(cfg["paths"]["reports"])
    reports.mkdir(parents=True, exist_ok=True)

    manifest = {
        "stage": stage,
        "project": cfg.get("project", {}).get("name"),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "config_path": cfg.get("_config_path"),
        "config_sha256": file_sha256(cfg.get("_config_path", "")),
        "inputs": inputs,
        "outputs": outputs,
    }
    if extra:
        manifest["details"] = extra

    out = reports / f"manifest_{stage}.json"
    out.write_text(json.dumps(manifest, indent=2))
    print(f"[manifest] {out}")
    return out


def banner(title: str) -> None:
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)
