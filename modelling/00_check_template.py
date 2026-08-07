#!/usr/bin/env python3
"""
00_check_template.py — verify the scaffold template before anything else.

Confirms the {} placeholder exists, reports where it sits in the sequence, and
checks that position against region.range in the config. A mismatch here
silently shifts every residue number downstream (SASA, RMSF, AIRs), so this is
worth.

Usage:
    python3 00_check_template.py config/hadv_c5_hvr7.yaml
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_config, load_template, banner

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else "config/hadv_c5_hvr7.yaml")
tmpl = load_template(cfg)

before = tmpl.split("{}")[0]
after = tmpl.split("{}")[1]
start = len(before) + 1                       # 1-based residue of graft start
query = cfg["region"]["query"]
rng = cfg["region"]["range"]

banner(f"Template check — {cfg['project']['name']}")
print(f"  scaffold length (excl. graft) : {len(before) + len(after)}")
print(f"  graft starts at residue       : {start}")
print(f"  config region.range           : {rng}")
print(f"  wild-type region query        : {len(query)} aa")
print(f"  full length with WT region    : {len(before) + len(query) + len(after)}")
print(f"  config protomer_length        : {cfg['target']['protomer_length']}")

ok = True
if start != rng[0]:
    print(f"\n  [MISMATCH] graft begins at {start} but region.range starts at {rng[0]}")
    ok = False
expected_len = rng[1] - rng[0] + 1
if len(query) != expected_len:
    print(f"  [NOTE] query is {len(query)} aa; region.range spans {expected_len}")
    print(f"         Fine if intentional, but confirm the graft replaces what you think.")

full = len(before) + len(query) + len(after)
if full != cfg["target"]["protomer_length"]:
    print(f"\n  [MISMATCH] template+query = {full} but protomer_length = "
          f"{cfg['target']['protomer_length']}")
    ok = False

print("\n  VERDICT:", "PASS" if ok else "FIX CONFIG BEFORE PROCEEDING")
sys.exit(0 if ok else 1)
