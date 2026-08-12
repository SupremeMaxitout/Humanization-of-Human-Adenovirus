# Humanisation of Human Adenovirus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21910027.svg)](https://doi.org/10.5281/zenodo.21910027)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![GROMACS 2023.2](https://img.shields.io/badge/GROMACS-2023.2-orange.svg)](https://www.gromacs.org/)
[![HADDOCK 2.5](https://img.shields.io/badge/HADDOCK-2.5-green.svg)](https://www.bonvinlab.org/)


**A reproducible computational pipeline for immune-evasive capsid engineering.**

Grafts human-derived sequence into the hypervariable regions (HVRs) of the
adenovirus hexon, then screens the resulting variants for structural quality,
immune recognition, and stability under physiological conditions.

Developed on HAdV-C5 hexon HVR7. Configuration-driven, so it retargets to other
adenovirus serotypes and other HVRs without editing code.

---

## Overview

```
   sequence                structure               interaction              dynamics
  ──────────              ───────────             ─────────────            ──────────
  proteome      ──▶  humanised trimer   ──▶   quality gate    ──▶   immune docking   ──▶   100 ns MD
  homology           AlphaFold3 /             geometry,             HADDOCK vs             stability
  search             ColabFold               confidence, SASA       component panel        screen

    STAGE 1              STAGE 1                  STAGE 2               STAGE 3            STAGE 4
```

| Stage | Directory | Purpose | Status |
|:-----:|-----------|---------|:------:|
| 1 | [`1_modelling/`](1_modelling/) | Proteome homology search → humanised trimer models → PDB | ✅ |
| 2 | [`2_validation/`](2_validation/) | Geometry, predictor confidence and SASA quality gate | ✅ |
| 3 | [`3_docking/`](3_docking/) | HADDOCK docking against an immune-component panel | ✅ |
| 4 | [`4_simulation/`](4_simulation/) | GROMACS 100 ns MD stability screen and analysis | ✅ |

Each stage runs standalone but shares one configuration file and one naming
convention. Directories are numbered so execution order is unambiguous.

---

## Quick start

```bash
git clone https://github.com/SupremeMaxitout/Humanization-of-Human-Adenovirus.git
cd Humanization-of-Human-Adenovirus
pip install -r requirements.txt
```

<details>
<summary><b>Stage 1 — Modelling</b></summary>

```bash
# human proteome (~110 MB uncompressed)
cd data
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz && cd ..

python3 1_modelling/00_check_template.py         config/hadv_c5_hvr7.yaml
python3 1_modelling/01_parse_proteome.py         config/hadv_c5_hvr7.yaml
python3 1_modelling/02_homology_search.py        config/hadv_c5_hvr7.yaml
python3 1_modelling/03_make_prediction_inputs.py config/hadv_c5_hvr7.yaml

# predict with AlphaFold3 server or ColabFold, then:
python3 1_modelling/04_unpack_predictions.py     config/hadv_c5_hvr7.yaml
python3 1_modelling/05_select_best_model.py      config/hadv_c5_hvr7.yaml
python3 1_modelling/06_cif_to_pdb.py             config/hadv_c5_hvr7.yaml
python3 1_modelling/07_make_docking_input.py     config/hadv_c5_hvr7.yaml
```
</details>

<details>
<summary><b>Stage 2 — Validation</b></summary>

```bash
python3 2_validation/01_geometry_check.py     config/hadv_c5_hvr7.yaml
python3 2_validation/02_confidence_metrics.py config/hadv_c5_hvr7.yaml
python3 2_validation/03_region_sasa.py        config/hadv_c5_hvr7.yaml
python3 2_validation/04_compile_report.py     config/hadv_c5_hvr7.yaml
```

Produces `reports/validation/validation_report.xlsx` and `passed_models.txt`,
the input list for stage 3.
</details>

<details>
<summary><b>Stage 3 — Docking</b></summary>

```bash
python3 3_docking/00_check_inputs.py config/hadv_c5_hvr7.yaml

TARGET_NAME=HD5 \
TARGET_PDB=$HOME/immune_system_components/HD5.pdb \
TARGET_AIR=$HOME/trimer_air/HD5_trimer_air.tbl \
bash 3_docking/dock_batch.sh

for f in $HOME/dock_HD5_run/jobs/*.pbs; do qsub "$f"; done

python3 3_docking/02_collect_results.py config/hadv_c5_hvr7.yaml
python3 3_docking/03_score_report.py    config/hadv_c5_hvr7.yaml
```

Repeat per immune component.
</details>

<details>
<summary><b>Stage 4 — Simulation</b></summary>

```bash
bash 4_simulation/generate_all_jobs.sh
bash 4_simulation/tools/resume_all.sh setup
bash 4_simulation/tools/check_status.sh
bash 4_simulation/tools/resume_all.sh prod

module load gromacs/2023.2-gpu
bash 4_simulation/tools/concat_trajectories.sh
bash 4_simulation/analyse_md.sh
```
</details>

---

## Documentation

| Document | Covers |
|---|---|
| [`docs/LESSONS.md`](docs/LESSONS.md) | **Read this first.** Every failure encountered building this, and its fix |
| [`docs/VALIDATION.md`](docs/VALIDATION.md) | What the quality gate checks, why, and how to set thresholds |
| [`docs/DOCKING.md`](docs/DOCKING.md) | Why HADDOCK, what the score means, why AIR files are the key assumption |
| [`docs/SIMULATION.md`](docs/SIMULATION.md) | MD parameters, troubleshooting, and how to read each output |

`LESSONS.md` documents roughly two weeks and 15,000 CPU-hours of failures — a
missing `refcoord-scaling = com`, a conjugate-gradient stage that could not
survive water constraints, and a "cleanup" step that manufactured the defects it
was meant to remove. Every one is now fixed in code or checked automatically.

---

## Retargeting

Copy `config/hadv_c5_hvr7.yaml`, change five things, run the same scripts:

```yaml
target:
  virus: HAdV-C2
  protomer_length: 939
  template: config/templates/hadv_c2_hexon_hvr1.txt   # scaffold with one {}

region:
  name: HVR1                                          # the region to humanise
  query: <wild-type sequence of that region>
  range: [136, 169]

regions:                                              # single source of truth
  HVR1: [136, 169]
  ...
```

`1_modelling/00_check_template.py` verifies the graft position matches the
residue numbering before you spend compute.

---

## Design decisions

**Multi-chain for MD, single-chain for docking.** The trimer is written twice.
HADDOCK treats each input as one rigid body; GROMACS needs three real chains.
Merging protomers for MD creates artificial peptide bonds across chain
junctions that stiffen the very interface an HVR flexibility study measures.

**No coordinate "cleanup" before simulation.** Model coordinates go into
`pdb2gmx` unmodified. Vacuum minimisation of a predicted structure manufactures
the defects it is meant to remove — measured here: 0 → 8 D-amino-acid centres
and 6 → 31 hard clashes on the same file.

**Region-restricted confidence, not just the global mean.** A model can score 90
pLDDT overall while the grafted loop is a low-confidence mess. The graft is the
subject of the study, so it is scored separately.

**Cluster-level docking scores.** The single best-scoring pose is largely noise.
Ranking uses the mean score of the best members of the best-populated cluster,
with cluster size reported as evidence.

**Thresholds fixed before results are seen.** The validation gate lives in the
config and is written into the report, so the applied criteria are recorded
alongside the numbers.

**Provenance on every step.** Each stage writes a manifest with timestamp, git
commit, config hash and input hashes.

---

## Limitations

Stated up front, because a screen that overclaims is worse than one that does not.

- **No experimental validation.** This prioritises candidates for wet-lab work;
  it does not predict immune evasion.
- **Single-replicate MD by default.** A screen, not converged sampling.
  Differences comparable to run-to-run noise are not real differences. Set
  `REPLICATES=3` for error bars.
- **Docking scores are not affinities.** HADDOCK ranks poses within one target;
  cross-target comparison is not supported, and HADDOCK 2.5 and 3 scores must
  never be pooled.
- **pLDDT is confidence, not accuracy.** A confidently predicted wrong structure
  scores well.
- **Rebuilt loops are models, not data.** Where a disordered HVR was rebuilt,
  the coordinates are predicted.

---

## Requirements

Python dependencies are in [`requirements.txt`](requirements.txt). External
software is not redistributed here.

| Software | Access |
|---|---|
| AlphaFold3 server | Free for academic use; **manual upload** — no public submission API |
| ColabFold | Open source, fully scriptable ([localcolabfold](https://github.com/YoshitakaMo/localcolabfold)) |
| HADDOCK 2.5 | On request from the Bonvin lab; licences are personal and non-transferable |
| HADDOCK3 | Open source, `pip install haddock3` |
| CNS | Free for non-profit use; register at [cns-online.org](http://cns-online.org/v1.3/) |
| GROMACS | Open source |

Stages 3 and 4 assume a PBS cluster (developed on NSCC ASPIRE 2A). SLURM
equivalents are documented in [`3_docking/README.md`](3_docking/README.md).

---

## Repository layout

```
config/            run configuration and scaffold templates
1_modelling/       candidate discovery and structure prediction
2_validation/      model quality gate
3_docking/         HADDOCK setup, execution and scoring
4_simulation/      GROMACS MD pipeline and analysis
docs/              methodology and lessons learned
data/              working directory (git-ignored)
reports/           spreadsheets, figures and provenance manifests
```

---

## Citation

If this pipeline contributes to your work, please cite this repository together
with the underlying tools: AlphaFold3, ColabFold, HADDOCK, CNS, GROMACS,
Biopython and gemmi.
