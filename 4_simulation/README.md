# 4 — Simulation

100 ns all-atom MD on each system under physiological conditions (310 K,
150 mM NaCl, explicit water), to test whether the engineered region destabilised
the trimer.

Input is the **multi-chain** structures from stage 1 — three real chains for a
homotrimer. The single-chain form is for docking only.

> Read [`docs/LESSONS.md`](../docs/LESSONS.md) before running anything. It
> documents roughly two weeks and 15,000 CPU-hours of failures, all now either
> fixed in code or checked automatically by the preflight. The single most
> important line in this stage is `refcoord-scaling = com`.

---

## Run

```bash
# 0. install
scp -r 4_simulation khan0023@asp2a-login-ntu01:~/md_pipeline
chmod +x ~/md_pipeline/*.sh ~/md_pipeline/tools/*.sh

# 1. list your systems
nano ~/md_pipeline/systems.conf

# 2. generate jobs
cd ~ && bash ~/md_pipeline/generate_all_jobs.sh

# 3. equilibrate (~1-2 h per system on one GPU)
bash ~/md_pipeline/tools/resume_all.sh setup

# 4. watch
bash ~/md_pipeline/tools/check_status.sh

# 5. production (2-3 submissions of 24 h each)
bash ~/md_pipeline/tools/resume_all.sh prod

# 6. join part files, then analyse
module load gromacs/2023.2-gpu
bash ~/md_pipeline/tools/concat_trajectories.sh
bash ~/md_pipeline/analyse_md.sh
```

## The three commands you will actually use

```bash
bash tools/check_status.sh          # where is everything?
bash tools/diagnose.sh <run>        # why did that one stop?
bash tools/resume_all.sh prod       # continue everything unfinished
```

`check_status.sh` distinguishes *queued*, *running* and *stopped*. The naive
`[OK]/[FAIL]` check it replaces reported healthy queued jobs as failures and
sent us chasing ghosts for days.

---

## Interrupted? Just resubmit.

Every stage is idempotent and every `mdrun` uses `-cpi`. Walltime kills, VPN
drops, cluster maintenance and node failures all resume from the last
checkpoint with no lost work.

```bash
bash tools/resume_all.sh setup   # or prod
```

Production reaching 100 ns in 2–3 submissions is **designed behaviour**, not
failure. 24 h is the ASPIRE2A maximum; 26 h is rejected. `-maxh 23.0` stops
gracefully and checkpoints.

Restart from scratch **only** if inputs or parameters changed — see
[`docs/SIMULATION.md`](../docs/SIMULATION.md) → "Resume vs restart".

---

## Your choices

### Force field

```bash
FORCE_FIELD=charmm27  bash generate_all_jobs.sh    # default
FORCE_FIELD=charmm36m bash generate_all_jobs.sh
```

`charmm27` ships with GROMACS and is what the reference runs used. `charmm36m`
is better for long disordered loops — which is what an HVR is — but is a
separate port from the MacKerell lab that you must install into `$GMXLIB`
yourself.

Whichever you choose: use the **same** one for every system in a comparison,
and state it accurately. They are different force fields, and labelling one as
the other is a real methods error.

### Replicates

```bash
REPLICATES=3 bash generate_all_jobs.sh
```

A single 100 ns run per system is a **screen**, not converged sampling.
Differences of the order of run-to-run noise are not real differences. Three
replicates with independent velocity seeds let you put error bars on the
per-region comparison — which is what a reviewer will ask for. It also triples
the compute.

Default is 1. Choose deliberately and say which you used.

---

## Files

```
4_simulation/
  systems.conf              which systems to run
  regions.conf              region boundaries, defined ONCE
  setup_md.sh               EM -> NVT -> warm-up -> NPT (resumable, preflight)
  generate_all_jobs.sh      builds PBS scripts from systems.conf
  job_setup.pbs.template    24 h, GPU, restart-safe
  job_prod.pbs.template     24 h, -maxh 23.0 -noappend -cpi
  analyse_md.sh             PBC-corrected analysis + figures
  mdp/                      em, nvt, npt_warmup, npt, prod
  tools/
    check_status.sh         stage reached, queue state, ns done
    diagnose.sh             why did this run stop?
    resume_all.sh           resubmit everything unfinished
    concat_trajectories.sh  join -noappend part files
```

---

## Output

```
~/md_runs/analysis/
  fig1_rmsd.png          did each system settle?
  fig2_rmsf.png          per-residue flexibility, regions shaded  <- KEY
  fig3_rg.png            trimer compactness
  fig4_region_bars.png   per-region RMSF vs reference  <- DECISION FIGURE
  summary.txt            RMSD mean, drift rate, Rg mean
```

A candidate passes the stability screen if RMSD plateaus near the reference,
per-region RMSF is not markedly higher, Rg stays flat, and the inter-protomer
H-bond count holds. Full interpretation in
[`docs/SIMULATION.md`](../docs/SIMULATION.md).

---

## Cost (per system, ~420k atoms, 1× A100)

| Phase | Time | Notes |
|---|---|---|
| Setup (EM + NVT + warm-up + NPT) | 1–2 h | was 10 h+ on CPU and hit walltime |
| Production 100 ns | 40–60 h | 2–3 × 24 h submissions, ~63 ns/day |

Six systems ≈ 19,000 SU.

---

## Cluster rules

- **Never run compute on a login node.** Three violations can block your
  account. `setup_md.sh` refuses to start on one.
- **Submit to `normal`, not `g1`.** PBS routes GPU jobs to `g1` itself;
  `qsub -q g1` is denied. `g1` also runs one job per user at a time, so jobs
  serialise — normal, not a fault.
- **`cd ~` before any `rm -rf`.** Deleting the directory you are standing in
  breaks `qsub` with `getcwd: No such file or directory`.
- **Never judge a run before `qstat` is empty.**
