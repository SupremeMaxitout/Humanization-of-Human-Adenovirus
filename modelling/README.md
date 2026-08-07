# Tool 1 — Candidate discovery and structure modelling

Human proteome → homologous candidate grafts → humanised trimer models →
PDB files ready for MD (multi-chain) and docking (single-chain).

Every stage is driven by one YAML config. Change the config, not the scripts.

---

## Pipeline

| Step | Script | In → Out |
|---|---|---|
| 00 | `00_check_template.py` | config → verify graft position vs residue numbering |
| 01 | `01_parse_proteome.py` | Ensembl FASTA → `ensembl_proteins.json` |
| 02 | `02_homology_search.py` | proteome + region query → leads `.xlsx` |
| 03 | `03_make_prediction_inputs.py` | leads → AF3 JSON **and** ColabFold FASTA |
| — | *(predict)* | AF3 server upload, or `colabfold_batch` |
| 04 | `04_unpack_predictions.py` | zips / CF output → uniform per-candidate dirs |
| 05 | `05_select_best_model.py` | seeds → best model + **all confidence metrics** |
| 06 | `06_cif_to_pdb.py` | mmCIF → PDB, 3 chains (**for MD**) |
| 07 | `07_make_docking_input.py` | 3 chains → 1 chain (**for HADDOCK**) |

```bash
python3 tool1_modelling/00_check_template.py        config/hadv_c5_hvr7.yaml
python3 tool1_modelling/01_parse_proteome.py        config/hadv_c5_hvr7.yaml
python3 tool1_modelling/02_homology_search.py       config/hadv_c5_hvr7.yaml
python3 tool1_modelling/03_make_prediction_inputs.py config/hadv_c5_hvr7.yaml
#   ... predict ...
python3 tool1_modelling/04_unpack_predictions.py    config/hadv_c5_hvr7.yaml
python3 tool1_modelling/05_select_best_model.py     config/hadv_c5_hvr7.yaml
python3 tool1_modelling/06_cif_to_pdb.py            config/hadv_c5_hvr7.yaml
python3 tool1_modelling/07_make_docking_input.py    config/hadv_c5_hvr7.yaml
```

Start with `--limit` on steps 01–02 to smoke-test in seconds.

---

## Multi-chain vs single-chain — read this once

Step 06 writes **3-chain** trimers. Step 07 writes **1-chain** trimers.

- **MD (GROMACS) uses the 3-chain files.**
- **Docking (HADDOCK) uses the 1-chain files.**

Feeding a single-chain file to `pdb2gmx` creates real peptide bonds across the
protomer junctions, artificially stiffening the inter-protomer interface — the
exact region an HVR flexibility study measures. That mistake produced three
unusable MD trajectories on this project before it was caught.

---

## Tuning the search

All in `config/*.yaml` under `homology:`

| Key | Effect |
|---|---|
| `window_size` | length of the candidate graft |
| `substitution_matrix` | BLOSUM45/62/80, PAM250 — 45 for remote homology, 80 for close |
| `alignment_mode` | `local` (default) or `global` |
| `score_threshold` | minimum alignment score |
| `charge_deviation` | tolerated \|Δ net charge\| vs the wild-type region |
| `match_polarity` | require the same charge sign |
| `ph` | pH for charge calculation |
| `dedupe_on` | collapse hits by `description`, `protein_ID`, or sequence |

**Why `local` is the default.** Globally aligning a 15-residue window against a
41-residue query charges the length difference as end gaps, so scores largely
track length rather than similarity. Local alignment compares like with like.
Set `global` to reproduce the original behaviour.

**Sanity check your thresholds.** `charge_deviation: 0.005` is extremely tight —
it demands near-identical net charge and will reject almost everything. Widen it
(0.5–1.0) if you get no hits, and say in your methods what you used.

---

## Graft length vs region length

If `window_size` (15) is shorter than the wild-type region (41 aa), the
humanised protomer is **shorter than wild type** — 900 vs 926 residues in the
worked example. That is a legitimate design choice, but it has consequences:

- residue numbering after the graft shifts between variants and wild type
- RMSF and SASA comparisons must be aligned on the scaffold, not raw numbering
- HADDOCK AIR definitions must use each model's own numbering

Step 03 prints the full-length range and warns when candidates differ. Decide
deliberately, and record it.

---

## Bugs fixed from the original scripts

| Original | Problem | Now |
|---|---|---|
| `script02` | `json.dump(..., indecnt=2)` — crashed; last record dropped | fixed |
| `script01` | `compare_with_QUERY` called twice per protein | once |
| `script01` | full alignment computed + printed per protein, then discarded | removed (dominant cost over 110k proteins) |
| `script01` | global alignment of 15-mer vs 41-mer | `local` default, configurable |
| `script05` | hardcoded `range(74)` | globs whatever is present |
| `script07` | output named `x.cif.pdb` | `.cif` stripped properly |
| `script08` | shifts hardcoded `{B:926, C:1852}` | derived from `protomer_length` |
| `script09` | HVR ranges duplicated and **disagreed** with the MD script | single source of truth in config |

That last one mattered: the SASA script used `HVR1 (136,169) / HVR7 (418,458)`
while the MD analysis used `HVR1 (137,181) / HVR7 (422,442)`. Two tools were
reporting on different residues under the same label.

---

## Provenance

Every step writes `reports/manifest_<stage>.json` with timestamp, git commit,
config SHA-256, input hashes, and counts. This is what makes the run
reproducible rather than merely scripted.
