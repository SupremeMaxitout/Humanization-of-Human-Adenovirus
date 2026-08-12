# 3 — Docking

Docks each validated model against an immune-component panel using HADDOCK,
and ranks the results by cluster score against the wild-type reference.

Input is the **single-chain** structures from stage 1 (`07_make_docking_input.py`)
filtered by `reports/validation/passed_models.txt` from stage 2.

```bash
python3 3_docking/00_check_inputs.py     config/hadv_c5_hvr7.yaml
python3 3_docking/01_generate_jobs.py    config/hadv_c5_hvr7.yaml
# submit the generated PBS jobs, wait, then:
python3 3_docking/02_collect_results.py  config/hadv_c5_hvr7.yaml
python3 3_docking/03_score_report.py     config/hadv_c5_hvr7.yaml
```

## Stages

| Script | Does |
|---|---|
| `00_check_inputs.py` | verifies HADDOCK, CNS, every target PDB + AIR file, and that ligands are single-chain |
| `01_generate_jobs.py` | writes one PBS job per ligand x target, for either engine |
| `02_collect_results.py` | rescues structures from `work_pbs/` into `results/` |
| `03_score_report.py` | cluster-level scoring, ranked vs wild type, to Excel |

## Config

The panel is fully generic — any number of targets, each with its own PDB and
AIR restraint file. Nothing about the immune components is hard-coded.

```yaml
paths:
  models_docking: data/models_docking    # SINGLE-CHAIN structures

docking:
  engine: haddock2.5                     # or haddock3
  haddock_dir: ~/software/haddock2.5-2025-08
  cns_exe: ~/software/cns_solve/bin/cns  # optional; verified if given
  base_dir: data/docking

  n_it0: 1000                            # rigid-body models
  n_it1: 200                             # semi-flexible refinement
  n_water: 200                           # final water-refined structures
  clust_cutoff: 7.5                      # FCC clustering
  clust_size: 4

  ncpu: 16
  mem: 64GB
  walltime: "08:00:00"

  targets:                               # one entry per immune component
    HD5:
      pdb: ~/immune_system_components/HD5.pdb
      air: ~/trimer_air/HD5_trimer_air.tbl
    FX:
      pdb: ~/immune_system_components/FX.pdb
      air: ~/trimer_air/FX_trimer_air.tbl
```

Run one target at a time with `--target HD5` on stages 01 and 02.

## Output

```
reports/docking/
  docking_report.xlsx
    Summary      ranked, all targets
    T_<target>   one sheet per target
    Full         every column
    HowToRead    what the numbers mean
  manifest_score.json
```

Ranking is by **cluster score** — the mean HADDOCK score of the best 4 members
of the best-scoring cluster — not the single best structure, which is mostly
noise. `cluster_size` and `cluster_sd` are reported alongside because a tight
cluster of 40 is far stronger evidence than a scattered cluster of 4.

`delta_vs_reference` is computed **within each target only**. HADDOCK scores are
not comparable across different targets, and scores from HADDOCK 2.5 and
HADDOCK3 must never be pooled.

See [`docs/DOCKING.md`](../docs/DOCKING.md) for why HADDOCK, what the score
means, and why AIR files are the biggest assumption in this stage.

## Requirements

Neither is redistributed here.

| Software | Access |
|---|---|
| HADDOCK 2.5 | on request from the Bonvin lab; licences are personal |
| HADDOCK3 | open source, `pip install haddock3` |
| CNS | free for non-profit, register at cns-online.org |
