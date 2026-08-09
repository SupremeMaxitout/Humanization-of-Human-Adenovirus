# Tool 2 — Validation

The quality gate between structure prediction (Tool 1) and docking (Tool 3).

Four stages, run in order. Each writes a CSV; the fourth merges them into one
spreadsheet and applies a pass/fail gate defined in the config.

```bash
python3 validation/01_geometry_check.py     config/hadv_c5_hvr7.yaml
python3 validation/02_confidence_metrics.py config/hadv_c5_hvr7.yaml
python3 validation/03_region_sasa.py        config/hadv_c5_hvr7.yaml
python3 validation/04_compile_report.py     config/hadv_c5_hvr7.yaml
```

Outputs land in `reports/validation/`:

| File | Contents |
|---|---|
| `geometry.csv` | D-centres, clashes, ring piercings, chain count |
| `confidence.csv` | pLDDT global and per region, pTM, ipTM, PAE |
| `sasa.csv` | per-region SASA and ratios to the reference |
| `validation_report.xlsx` | all of the above plus the verdict |
| `passed_models.txt` | model list that feeds Tool 3 |
| `manifest_*.json` | provenance for each stage |

## What each stage is for

**01 — geometry.** Catches the three defect classes that no minimiser can
repair: inverted stereocentres (D-amino acids), hard steric overlaps below
1.2 Å, and bonds threaded through aromatic rings. Each survives energy
minimisation and restrained equilibration, then destroys the run once dynamics
are unrestrained. Mild 1.8–2.3 Å polar contacts are reported but not failed —
they are normal in predicted structures.

**02 — confidence.** Global pLDDT plus pLDDT restricted to each region. The
region number is the one that matters: a model can score 90 overall while the
grafted loop is a low-confidence mess, and the graft is the entire subject of
the study. Also collects pTM, ipTM (trimer assembly confidence), PAE, and the
predictor's own clash flag. Works with AlphaFold3 and ColabFold output.

**03 — SASA.** A graft only evades an antibody if it is still presented on the
surface. Reported as a ratio to the wild-type reference; ~1.0 means
presentation is preserved.

**04 — report.** Merges everything, applies the gate, writes the spreadsheet
and the pass list.

## Options

```bash
--models DIR      model directory (default: paths.models in config)
--json-dir DIR    predictor confidence JSONs (stage 02)
--reference NAME  SASA baseline (default: validation.reference in config)
--n-points N      Shrake-Rupley sphere points, stage 03 (default 100)
```

## Config

```yaml
target:
  protomer_length: 952      # folds chain B/C numbering onto protomer 1
region:
  name: HVR7                # the region being humanised
regions:                    # single source of truth for boundaries
  HVR1: [136, 169]
  ...
paths:
  models: data/models
  predictions: data/predictions
validation:
  reference: wild_type
  thresholds:               # all optional; see docs/VALIDATION.md
    min_mean_plddt: 70.0
    min_region_plddt: 70.0
    min_iptm: 0.60
```

Thresholds are set in the config and applied uniformly, so the gate is fixed
before results are seen.

## Notes

Missing metrics are reported as blanks, not failures. If no confidence JSONs
are supplied, stage 02 still reports pLDDT from the B-factor column and the
gate simply skips the checks it cannot evaluate.

Residue numbering is folded onto protomer 1, so the same region boundaries work
on three-chain models and merged single-chain models alike.
