# Humanisation of Human Adenovirus

A reproducible computational pipeline for **immune-evasive capsid engineering**:
grafting human-derived sequence into the hypervariable regions (HVRs) of the
adenovirus hexon, then screening the resulting variants for structural quality,
immune recognition, and stability.

Developed on HAdV-C5 hexon HVR7, written to generalise across adenovirus
serotypes and any HVR via configuration.

---

## Pipeline

| | Stage | Does | Status |
|---|---|---|---|
| 1 | **[Modelling](1_modelling/)** | proteome homology search → humanised trimer models → PDB | ✅ |
| 2 | **[Validation](2_validation/)** | geometry + AlphaFold confidence + SASA gate | ✅ |
| 3 | **[Docking](3_docking/)** | HADDOCK against an immune-component panel | 🔨 |
| 4 | **[Simulation](4_simulation/)** | GROMACS MD stability screen + analysis | 🔨 |

Each stage stands alone but shares one config and one naming convention.
Folders are numbered so the execution order is unambiguous.

---

## Quick start

```bash
git clone https://github.com/SupremeMaxitout/Humanization-of-Human-Adenovirus.git
cd Humanization-of-Human-Adenovirus
pip install -r requirements.txt

# get the human proteome (~110 MB uncompressed)
cd data
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
cd ..
```

**Stage 1 — modelling**

```bash
python3 1_modelling/00_check_template.py       config/hadv_c5_hvr7.yaml
python3 1_modelling/01_parse_proteome.py       config/hadv_c5_hvr7.yaml
python3 1_modelling/02_homology_search.py      config/hadv_c5_hvr7.yaml
python3 1_modelling/03_make_prediction_inputs.py config/hadv_c5_hvr7.yaml
# → predict with AlphaFold3 server or ColabFold, then:
python3 1_modelling/04_unpack_predictions.py   config/hadv_c5_hvr7.yaml
python3 1_modelling/05_select_best_model.py    config/hadv_c5_hvr7.yaml
python3 1_modelling/06_cif_to_pdb.py           config/hadv_c5_hvr7.yaml
python3 1_modelling/07_make_docking_input.py   config/hadv_c5_hvr7.yaml
```

**Stage 2 — validation**

```bash
python3 2_validation/01_geometry_check.py      config/hadv_c5_hvr7.yaml
python3 2_validation/02_confidence_metrics.py  config/hadv_c5_hvr7.yaml
python3 2_validation/03_region_sasa.py         config/hadv_c5_hvr7.yaml
python3 2_validation/04_compile_report.py      config/hadv_c5_hvr7.yaml
```

Produces `reports/validation/validation_report.xlsx` and `passed_models.txt`,
the input list for stage 3.

---

## Retargeting to another serotype or region

Copy `config/hadv_c5_hvr7.yaml`, change five things, run the same scripts:

```yaml
target:
  virus: HAdV-C2
  protomer_length: 939
  template: config/templates/hadv_c2_hexon_hvr1.txt   # scaffold with one {}
region:
  name: HVR1
  query: <wild-type sequence of that region>
  range: [136, 169]
regions:            # all region boundaries — single source of truth
  HVR1: [136, 169]
  ...
```

`1_modelling/00_check_template.py` verifies the graft position matches the
residue numbering before you spend compute.

---

## Design decisions worth knowing

**Multi-chain for MD, single-chain for docking.** The trimer is written twice.
HADDOCK wants one rigid body; GROMACS needs three real chains. Merging them for
MD creates artificial peptide bonds across protomer junctions that stiffen the
very interface an HVR flexibility study measures.

**No coordinate "cleanup" before simulation.** Model coordinates go into
`pdb2gmx` unmodified. Vacuum minimisation of a predicted structure manufactures
the defects it is meant to remove — measured on this project: 0 → 8
D-amino-acid centres and 6 → 31 hard clashes on the same file.

**Region-restricted confidence, not just the global mean.** A model can score 90
pLDDT overall while the grafted loop specifically is a low-confidence mess. The
graft is the entire subject of the study, so it is scored separately.

**One source of truth for region boundaries.** Defined once in the config.
Previously duplicated across scripts that quietly disagreed by several residues.

**Thresholds fixed before results are seen.** The validation gate lives in the
config and is written into the report, so the applied criteria are recorded
alongside the numbers.

**Provenance on every step.** Each stage writes a manifest with timestamp, git
commit, config hash, and input hashes.

---

## Limitations — stated up front

- **No experimental validation.** This prioritises candidates for wet-lab work;
  it does not predict immune evasion.
- **Single-replicate MD.** A screen, not converged sampling. Differences on the
  order of run-to-run noise are not real differences.
- **Docking scores are not affinities.** HADDOCK scores rank poses within a
  target; cross-target comparison is soft.
- **pLDDT is confidence, not accuracy.** A confidently predicted wrong structure
  scores well.
- **Rebuilt loops are models, not data.** Where a disordered HVR was rebuilt,
  the coordinates are predicted.

---

## Third-party requirements

| Software | Access |
|---|---|
| AlphaFold3 server | free, academic, **manual upload** — no public submission API |
| ColabFold | open source, fully scriptable ([localcolabfold](https://github.com/YoshitakaMo/localcolabfold)) |
| HADDOCK3 | open source; **requires CNS**, free for non-profit, registered separately |
| HADDOCK2.5 | on request from the Bonvin lab; licences are personal and not redistributable |
| GROMACS | open source |

No licensed software, credentials, or model weights are redistributed here.

---

## Repository layout

```
config/            run configs + scaffold templates
1_modelling/       candidate discovery and structure prediction
2_validation/      model quality gate
3_docking/         HADDOCK setup and scoring      (in progress)
4_simulation/      GROMACS MD pipeline            (in progress)
docs/              methodology and lessons learned
data/              working directory (git-ignored)
reports/           spreadsheets, figures, manifests
```

---

## Documentation

| Doc | Covers |
|---|---|
| [docs/VALIDATION.md](docs/VALIDATION.md) | what is checked, why, and how to choose thresholds |

---

## Citation

If this pipeline contributes to your work, please cite this repository and the
underlying tools: AlphaFold3, ColabFold, HADDOCK, GROMACS, Biopython, gemmi.
