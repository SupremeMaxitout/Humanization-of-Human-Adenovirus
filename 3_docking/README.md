# 3 — Docking

Docks each validated model against an immune-system component using HADDOCK,
then ranks the results by cluster score against the wild-type reference.

Input is the **single-chain** structures from stage 1
(`07_make_docking_input.py`), optionally filtered by
`reports/validation/passed_models.txt` from stage 2.

---

## Which script to use

`dock_batch.sh` is the **reference implementation** — the bash generator that
actually ran this project's docking campaign on NSCC ASPIRE 2A, generalised so
any target can be supplied. It has no dependencies beyond bash, so it runs on a
login node as-is. Use it for HADDOCK 2.5.

The Python scripts handle the parts bash is poor at: pre-flight validation,
rescuing results, and parsing hundreds of PDB headers into a ranked spreadsheet.
They run once, locally or on a login node, not on compute.

| Script | Language | Runs where | Purpose |
|---|---|---|---|
| `00_check_inputs.py` | Python | login node | verify HADDOCK, CNS, targets, AIRs, single-chain ligands |
| **`dock_batch.sh`** | **bash** | **login node → PBS** | **generate + submit docking jobs (HADDOCK 2.5)** |
| `01_generate_jobs.py` | Python | login node | alternative generator; also the HADDOCK3 path |
| `02_collect_results.py` | Python | login node | rescue structures from `work_pbs/` |
| `03_score_report.py` | Python | anywhere | cluster scoring → Excel |

The docking itself always runs as PBS batch jobs on compute nodes. Both
generators emit the same kind of PBS script; only the templating language
differs.

---

## Typical run

```bash
# 1. pre-flight
python3 3_docking/00_check_inputs.py config/hadv_c5_hvr7.yaml

# 2. generate jobs — ONE target per invocation
TARGET_NAME=HD5 \
TARGET_PDB=$HOME/immune_system_components/HD5.pdb \
TARGET_AIR=$HOME/trimer_air/HD5_trimer_air.tbl \
LIGAND_DIR=$HOME/trimers \
BASE_DIR=$HOME/dock_HD5_run \
bash 3_docking/dock_batch.sh

# 3. submit
for f in $HOME/dock_HD5_run/jobs/*.pbs; do qsub "$f"; done
qstat -u $USER

# 4. collect anything the jobs missed
python3 3_docking/02_collect_results.py config/hadv_c5_hvr7.yaml

# 5. score
python3 3_docking/03_score_report.py config/hadv_c5_hvr7.yaml
```

Repeat step 2–3 for each immune component: HD5, HD6, MARCO, Factor X,
neutralising antibody, or whatever panel your system needs.

---

## Configuring `dock_batch.sh`

Everything is an environment variable, so the script stays dependency-free on a
login node (no python, no `yq` needed to read YAML).

| Variable | Default | Meaning |
|---|---|---|
| `TARGET_NAME` | **required** | short label, e.g. `HD5` |
| `TARGET_PDB` | **required** | receptor structure |
| `TARGET_AIR` | **required** | AIR restraint file for this target |
| `LIGAND_DIR` | `$HOME/trimers` | single-chain models to dock |
| `BASE_DIR` | `$HOME/dock_<TARGET>_run` | output root |
| `HADDOCK_DIR` | `$HOME/software/haddock2.5-2025-08` | HADDOCK install |
| `PASSED_LIST` | unset | restrict to models listed in this file |
| `N_IT0` / `N_IT1` / `N_WATER` | `1000` / `200` / `200` | sampling |
| `CLUST_CUTOFF` / `CLUST_SIZE` | `7.5` / `4` | FCC clustering |
| `NCPU` / `MEM` / `WALLTIME` / `QUEUE` | `16` / `64GB` / `08:00:00` / `normal` | resources |

To dock only models that passed validation:

```bash
PASSED_LIST=reports/validation/passed_models.txt \
TARGET_NAME=HD5 ... bash 3_docking/dock_batch.sh
```

**On sampling:** the defaults are HADDOCK's recommended production values. The
pilot runs for this project used `N_IT0=50 N_IT1=25 N_WATER=25`, which is far
cheaper but gives too few models for cluster statistics to mean much. Choose
deliberately and state which you used.

The Python stages read the `docking:` block of the YAML config — see
`config/hadv_c5_hvr7.yaml` and the example in `docs/DOCKING.md`.

---

## Other schedulers

`dock_batch.sh` emits PBS (Torque/OpenPBS), because that is what ASPIRE 2A runs.
On a SLURM cluster, replace the header block in the generated jobs:

| PBS | SLURM |
|---|---|
| `#PBS -N name` | `#SBATCH --job-name=name` |
| `#PBS -q normal` | `#SBATCH --partition=normal` |
| `#PBS -l select=1:ncpus=16:mem=64GB` | `#SBATCH --nodes=1 --cpus-per-task=16 --mem=64G` |
| `#PBS -l walltime=08:00:00` | `#SBATCH --time=08:00:00` |
| `#PBS -j oe` / `#PBS -o path` | `#SBATCH --output=path` |
| `qsub job.pbs` | `sbatch job.sh` |
| `qstat -u $USER` | `squeue -u $USER` |

Everything below the header — the HADDOCK invocation, the fixes, the collection
logic — is scheduler-independent.

---

## Output

```
$BASE_DIR/
  jobs/                 generated PBS scripts
  logs/                 per-ligand job output
  work_pbs/run_<lig>/   HADDOCK working directories
  results/<lig>/
    it0/                rigid-body structures
    it1/                water-refined structures  ← the final models
    analysis/           cluster output

reports/docking/
  docking_report.xlsx
    Summary      ranked, all targets
    T_<target>   one sheet per target
    Full         every column
    HowToRead    what the numbers mean
```

Ranking is by **cluster score** — the mean HADDOCK score of the best 4 members
of the best-scoring cluster — not the single best structure, which is mostly
noise. `cluster_size` and `cluster_sd` sit alongside it, because a tight cluster
of 40 is far stronger evidence than a scattered cluster of 4.

`delta_vs_reference` is computed **within each target only**. HADDOCK scores are
not comparable across targets, and 2.5 and HADDOCK3 scores must never be pooled.

See [`docs/DOCKING.md`](../docs/DOCKING.md) for why HADDOCK, what the score
means, and why AIR files are the largest assumption in this stage.

---

## Requirements

Neither is redistributed here.

| Software | Access |
|---|---|
| HADDOCK 2.5 | on request from the Bonvin lab; licences are personal |
| HADDOCK3 | open source, `pip install haddock3` |
| CNS | free for non-profit, register at <http://cns-online.org/v1.3/> |
