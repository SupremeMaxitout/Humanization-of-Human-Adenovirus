"""
common.py — shared helpers for Tool 2 (validation).

Mirrors the conventions used in tool 1: one YAML config is the single source of
truth for region boundaries, thresholds and paths, and every stage writes a
manifest so a result can be traced back to the inputs that produced it.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML missing.  pip install -r requirements.txt")


# --------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------

# Defaults applied when the config does not override them. Every one of these
# is a decision, not a fact — see docs/VALIDATION.md for the reasoning.
DEFAULT_THRESHOLDS = {
    "min_mean_plddt": 70.0,       # whole model; <70 is "low confidence" per AF
    "min_region_plddt": 70.0,     # the grafted region specifically
    "min_iptm": 0.60,             # inter-chain confidence for the trimer
    "min_ptm": 0.50,
    "max_interface_pae": 10.0,    # Angstrom, mean PAE across chain pairs
    "max_d_centres": 0,
    "max_severe_clashes": 0,      # heavy-atom pairs < 1.2 A
    "max_ring_piercings": 0,
    "expected_chains": 3,
    "max_fraction_disordered": 0.30,
}


def load_config(path: str | Path) -> dict:
    """Load the run config and fill in validation defaults."""
    path = Path(path)
    if not path.is_file():
        sys.exit(f"[FATAL] config not found: {path}")
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}

    thr = dict(DEFAULT_THRESHOLDS)
    thr.update((cfg.get("validation") or {}).get("thresholds") or {})
    cfg.setdefault("validation", {})["thresholds"] = thr
    return cfg


def get_regions(cfg: dict) -> dict:
    """
    Region boundaries, e.g. {"HVR1": (137, 181), ...}.

    Single source of truth: defined once in the config so scripts cannot
    silently disagree about where a region starts.
    """
    raw = cfg.get("regions") or {}
    if not raw:
        sys.exit("[FATAL] config has no 'regions:' block — cannot score by region.")
    return {k: (int(v[0]), int(v[1])) for k, v in raw.items()}


def get_target_region(cfg: dict) -> str | None:
    """The region actually being humanised (config: region.name)."""
    return (cfg.get("region") or {}).get("name")


def protomer_length(cfg: dict) -> int | None:
    val = (cfg.get("target") or {}).get("protomer_length")
    return int(val) if val else None


# --------------------------------------------------------------------------
# residue numbering across a homotrimer
# --------------------------------------------------------------------------

def normalise_resnum(resnum: int, plen: int | None) -> int:
    """
    Map a residue number onto protomer-1 numbering.

    Region boundaries are defined for ONE protomer. Depending on how a model was
    written, chains B and C may be numbered 1..N again (separate chains) or
    offset by N and 2N (merged single chain). This collapses both onto 1..N so
    a region lookup works either way.
    """
    if not plen or plen <= 0:
        return resnum
    r = ((resnum - 1) % plen) + 1
    return r


def in_region(resnum: int, bounds: tuple[int, int], plen: int | None) -> bool:
    lo, hi = bounds
    return lo <= normalise_resnum(resnum, plen) <= hi


# --------------------------------------------------------------------------
# provenance
# --------------------------------------------------------------------------

def sha256(path: str | Path, blocks: int = 1 << 16) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(blocks), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "unknown"


def write_manifest(out_dir: str | Path, stage: str, cfg_path: str | Path,
                   inputs: list, extra: dict | None = None) -> Path:
    """Record what produced a result, so it can be traced later."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    man = {
        "stage": stage,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "config": str(cfg_path),
        "config_sha256": sha256(cfg_path) if Path(cfg_path).is_file() else None,
        "n_inputs": len(inputs),
        "inputs": [str(p) for p in inputs[:200]],
    }
    if extra:
        man.update(extra)
    dest = out_dir / f"manifest_{stage}.json"
    with open(dest, "w") as fh:
        json.dump(man, fh, indent=2)
    return dest


# --------------------------------------------------------------------------
# misc
# --------------------------------------------------------------------------

def model_id(path: str | Path) -> str:
    """
    Stable short name for a model file.

    Strips the extension chains that AF3 and ColabFold leave behind, so
    'fold_013_ensp...__model_0.cif.pdb' and its ColabFold equivalent both
    reduce to something usable as a spreadsheet key.
    """
    name = Path(path).name
    for suf in (".cif.pdb", ".pdb", ".cif"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    for marker in ("_unrelaxed", "_relaxed", "_single_chain"):
        if marker in name:
            name = name.split(marker)[0]
    return name


def banner(text: str) -> None:
    print("=" * 70)
    print(f"  {text}")
    print("=" * 70)
