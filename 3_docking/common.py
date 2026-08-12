"""
common.py — shared helpers for Tool 3 (docking).

Same conventions as Tools 1 and 2: one YAML config is the single source of
truth, and every stage writes a manifest so a result traces back to the inputs
and code that produced it.

Engine support
--------------
Both HADDOCK 2.5 and HADDOCK3 are supported, and every output records which
engine produced it (`engine` column, `engine` field in manifests, and the
`Engine` sheet in the report). They are NOT interchangeable:

  HADDOCK 2.5  run.param + run.cns, CNS-driven, the classic interface.
               This is the engine the reference project used.
  HADDOCK3     modular TOML workflow (`haddock3 run.toml`), different module
               names and output layout.

Scores from the two engines are not directly comparable. Rank within one engine
and say which one you used.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML missing.  pip install -r requirements.txt")


ENGINES = ("haddock2.5", "haddock3")

DEFAULT_DOCKING = {
    "engine": "haddock2.5",
    "n_it0": 1000,        # rigid-body models generated
    "n_it1": 200,         # semi-flexible refinement
    "n_water": 200,       # water refinement — these are the final structures
    "ncpu": 16,
    "mem": "64GB",
    "walltime": "08:00:00",
    "clust_cutoff": 7.5,  # FCC cutoff for clustering
    "clust_size": 4,      # minimum members for a cluster to count
}


def load_config(path: str | Path) -> dict:
    path = Path(path)
    if not path.is_file():
        sys.exit(f"[FATAL] config not found: {path}")
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}

    dock = dict(DEFAULT_DOCKING)
    dock.update(cfg.get("docking") or {})
    if dock["engine"] not in ENGINES:
        sys.exit(f"[FATAL] docking.engine must be one of {ENGINES}, "
                 f"got '{dock['engine']}'")
    cfg["docking"] = dock
    return cfg


def get_targets(cfg: dict) -> dict:
    """
    The immune-component panel.

    Generic by design: any number of targets, each a PDB plus its own AIR
    restraint file. Nothing about the panel is hard-coded — swap in whatever
    receptor set your system needs.

        docking:
          targets:
            HD5:
              pdb: immune_system_components/HD5.pdb
              air: trimer_air/HD5_trimer_air.tbl
            FX:
              pdb: immune_system_components/FX.pdb
              air: trimer_air/FX_trimer_air.tbl
    """
    targets = (cfg.get("docking") or {}).get("targets") or {}
    if not targets:
        sys.exit("[FATAL] config has no docking.targets — nothing to dock against.\n"
                 "        Define at least one target with 'pdb:' and 'air:'.")
    out = {}
    for name, spec in targets.items():
        if not isinstance(spec, dict) or "pdb" not in spec or "air" not in spec:
            sys.exit(f"[FATAL] target '{name}' needs both 'pdb:' and 'air:'")
        out[name] = {"pdb": Path(spec["pdb"]).expanduser(),
                     "air": Path(spec["air"]).expanduser()}
    return out


def get_ligand_dir(cfg: dict) -> Path:
    """
    Directory of structures to dock (the humanised models).

    Docking uses the SINGLE-CHAIN form written by Tool 1 stage 07: HADDOCK
    treats each input as one rigid body, so a three-chain trimer would be
    docked as three separate molecules. The multi-chain form is for MD only.
    """
    paths = cfg.get("paths") or {}
    return Path(paths.get("models_docking", "data/models_docking")).expanduser()


# --------------------------------------------------------------------------
# HADDOCK output parsing
# --------------------------------------------------------------------------

# HADDOCK 2.5 file.list lines look like:
#   "PREVIT:complex_1w.pdb"  { -142.318 }
FILE_LIST_RE = re.compile(r'"([^"]+)".*?\{\s*(-?[\d.]+)\s*\}')


def parse_file_list(path: Path) -> list[tuple[str, float]]:
    """Parse a HADDOCK file.list into [(structure, score), ...], best first."""
    out = []
    try:
        text = Path(path).read_text()
    except OSError:
        return out
    for line in text.splitlines():
        m = FILE_LIST_RE.search(line)
        if m:
            out.append((m.group(1).split(":")[-1], float(m.group(2))))
    return out


def parse_haddock_pdb_energies(path: Path) -> dict:
    """
    Pull energy terms from the REMARK header of a HADDOCK output PDB.

    HADDOCK writes lines such as:
        REMARK energies: <total>, <bonds>, <angles>, ..., <vdw>, <elec>, ...
        REMARK Desolvation energy: <value>
        REMARK buried surface area: <value>
    Field positions vary between versions, so this reads defensively and simply
    omits anything it cannot find rather than guessing.
    """
    vals: dict[str, float] = {}
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("REMARK"):
                    if line.startswith(("ATOM", "HETATM")):
                        break          # header is over
                    continue
                low = line.lower()
                if "energies:" in low:
                    nums = re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?",
                                      line.split(":", 1)[1])
                    if len(nums) >= 7:
                        vals["Evdw"] = float(nums[5])
                        vals["Eelec"] = float(nums[6])
                elif "desolvation energy" in low:
                    nums = re.findall(r"-?\d+\.?\d*", line)
                    if nums:
                        vals["Edesolv"] = float(nums[-1])
                elif "buried surface area" in low:
                    nums = re.findall(r"-?\d+\.?\d*", line)
                    if nums:
                        vals["BSA"] = float(nums[-1])
                elif "air energy" in low:
                    nums = re.findall(r"-?\d+\.?\d*", line)
                    if nums:
                        vals["Eair"] = float(nums[-1])
    except OSError:
        pass
    return vals


def haddock_score(vals: dict) -> float | None:
    """
    Standard HADDOCK score for water-refined models:

        1.0*Evdw + 0.2*Eelec + 1.0*Edesolv + 0.1*Eair

    Used only when a score is not already available from file.list. See
    docs/DOCKING.md for what the weights mean and what the score is not.
    """
    if "Evdw" not in vals or "Eelec" not in vals:
        return None
    return round(1.0 * vals["Evdw"]
                 + 0.2 * vals["Eelec"]
                 + 1.0 * vals.get("Edesolv", 0.0)
                 + 0.1 * vals.get("Eair", 0.0), 3)


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


def write_manifest(out_dir, stage: str, cfg_path, inputs: list,
                   extra: dict | None = None) -> Path:
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


def model_id(path) -> str:
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
