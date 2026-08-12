# Simulation methodology

Covers parameters, troubleshooting and analysis for stage 4. Read
[`LESSONS.md`](LESSONS.md) first — it explains *why* several of these
parameters are non-negotiable.

---

## The five-stage protocol

| Stage | File | Length | Restrained | Purpose |
|---|---|---|---|---|
| EM | `em.mdp` | to convergence | n/a | remove steric strain from the model |
| NVT | `nvt.mdp` | 100 ps | yes | bring to 310 K at fixed volume |
| NPT warm-up | `npt_warmup.mdp` | 50 ps @ 0.5 fs | yes | let the barostat act *gently* the first time |
| NPT | `npt.mdp` | 1 ns @ 1 fs | yes | settle density and pressure |
| Production | `prod.mdp` | 100 ns @ 2 fs | **no** | the only trajectory you analyse |

The warm-up exists because the first barostat contact is where model-derived
structures fail. Removing it saves ~30 min and risks a whole 24 h job.

---

## Parameters you will legitimately change

**`nsteps`** — simulation length. `nsteps × dt = total time`;
50,000,000 × 0.002 ps = 100 ns.
- 100 ns: adequate for a relative stability screen.
- 200–500 ns: needed for converged loop sampling or claims about a specific
  loop conformation.
- 3 × 100 ns replicates beats 1 × 300 ns for statistical confidence. Use
  `REPLICATES=3`.

**`ref-t`** — 310 K (physiological). Use 298 K only to match a
room-temperature experiment.

**Salt** (`SALT_CONC`) — 0.15 M NaCl approximates extracellular fluid, correct
for a capsid protein exposed to blood. Use 0.10 M for intracellular contexts.

**`BOX_PADDING`** — 1.5 nm minimum. Below ~1.2 nm the trimer interacts with its
own periodic image. Increase to 2.0 nm if loops extend far, but cost scales with
box volume.

**`dt`** — 2 fs with `constraints = h-bonds` for production; 1 fs for restrained
equilibration; 0.5 fs only in the warm-up. Never 2 fs without h-bond
constraints.

**`nstxout-compressed`** — 25,000 steps = 50 ps/frame = 2000 frames over 100 ns.

**Force field** — `charmm27` (ships with GROMACS, what the reference runs used)
or `charmm36m` (better for long disordered loops; separate install into
`$GMXLIB`). They are **different force fields**. Use one consistently across a
comparison and name it accurately.

---

## Parameters you should NOT change

| Parameter | Value | Why |
|---|---|---|
| `refcoord-scaling` | `com` | Mandatory with `-DPOSRES` + barostat. LESSONS #1. |
| `pcoupl` | `C-rescale` | Berendsen does not sample the correct NPT ensemble; Parrinello-Rahman is unstable far from equilibrium |
| `tcoupl` | `V-rescale` | Correct canonical ensemble; Berendsen has the flying-ice-cube artifact |
| `integrator` (EM) | `steep` | `cg` cannot handle water constraint failures. LESSONS #2. |
| `cutoff-scheme` | `Verlet` | the only supported scheme in modern GROMACS |
| `coulombtype` | `PME` | cutoff electrostatics are wrong for charged proteins |
| `rcoulomb` / `rvdw` | 1.2 nm | matched to CHARMM parameterisation |

---

## Troubleshooting

### Triage order

```bash
bash tools/check_status.sh        # where is everything?
bash tools/diagnose.sh <run>      # why did THIS one stop?
```

1. **Is it still running?** `qstat -u $USER`. A missing `.gro` on a queued job
   means nothing.
2. **Killed, or crashed?** `grep -i walltime *.log`. A walltime kill is not a
   physics problem.
3. **Did a tool refuse to start?** `grep -i "Fatal error" grompp_*.log`. Missing
   files and bad mdp options live here, and they are the most common cause.
4. **Only then** look at the structure.

### Where each error lives

| File | Contains |
|---|---|
| `setup_pbs.log` / `prod_pbs.log` | PBS stdout; walltime kills appear here |
| `pdb2gmx.log` | residue/atom recognition, chain detection |
| `grompp_*.log` | **mdp errors, missing files, topology mismatches** — check first |
| `mdrun_*.log` | LINCS warnings, water settle failures, segfaults |
| `em.log` | minimisation convergence, final Fmax |
| `npt*.log` / `prod*.log` | step-by-step energies; last `Step Time` block is progress |
| `step*.pdb` | coordinates dumped just before a crash — **check timestamps, these go stale** |

### Symptom → fix

| Symptom | Cause | Fix |
|---|---|---|
| `PBS: job killed: walltime` | ran out of clock | resubmit; `-cpi` resumes |
| `Run time exceeded ... will terminate` | `-maxh` graceful stop | resubmit |
| LINCS explosion at ~800–1500 steps of restrained NPT | missing `refcoord-scaling = com` | add it, delete `npt.tpr`, resubmit |
| `Epot = inf`, `Fmax = 2e+06` in EM | overlapping atoms in the **input** | fix the structure; do not minimise harder, never relax in vacuum |
| `One or more water molecules can not be settled` | local bad contact | `diagnose.sh` section 5 maps atom → residue |
| `Minimizer 'cg' can not handle constraint failures` | a cg stage was re-introduced | remove it |
| `Cannot do appending ... checksum` | PBS overwrote GROMACS's log | `mv prod.log prod_old.log`, use `-noappend` |
| `File 'npt_warmup.mdp' does not exist` | file never reached the cluster | check all five mdp files are present |
| `Access to queue is denied` | submitted with `-q g1` | use `-q normal` |
| `getcwd: No such file or directory` | current directory was deleted | `cd ~` |
| Segfault with no message | physics blow-up hidden by GPU offload | re-run with `USE_GPU=0` for readable LINCS output |

### Resume vs restart

**Resume** (the default — never loses completed work): walltime kills, VPN
drops, node failures, cluster maintenance.

```bash
bash tools/resume_all.sh setup    # or prod
```

**Restart from scratch** only when inputs or parameters changed: you edited an
mdp (delete the affected `.tpr`), changed the input PDB, force field, box size
or salt.

```bash
cd ~                      # never delete the directory you are standing in
rm -rf ~/md_runs/<run>
bash generate_all_jobs.sh
qsub ~/md_runs/<run>/job_setup.pbs
```

**Never** re-run `grompp` for production between restarts — a new `prod.tpr`
breaks checkpoint continuity.

---

## Analysis

```bash
bash tools/concat_trajectories.sh    # join -noappend part files FIRST
bash analyse_md.sh
SKIP_NS=20 bash analyse_md.sh        # discard 20 ns of equilibration (default)
```

### The raw files

| File | Use |
|---|---|
| `prod*.xtc` | coordinates, 50 ps/frame — every structural metric |
| `prod_full.xtc` | the joined trajectory; **this is what analysis reads** |
| `prod*.edr` | energies, T, P, volume — verify the thermostat behaved |
| `prod*.cpt` | checkpoint — never delete mid-run |
| `prod*.gro` | final frame; only written on genuine completion |
| `prod.tpr` | binary run input, required by every analysis command |

### The metrics

**RMSD** (`fig1_rmsd.png`) — did the fold hold? Rises then plateaus within
10–20 ns is normal; 1–3 Å is a stable globular protein. Still climbing at 50 ns
means not converged or genuinely unfolding. A sudden vertical jump means PBC
correction failed, not biology.

**RMSF** (`fig2_rmsf.png`) — **the key plot.** Per-residue fluctuation after
discarding equilibration. Compare each engineered region against the reference
*at the same residue positions*. Similar or lower = structurally tolerated.
Markedly higher (roughly >1.5×) = the graft is destabilising. A change *outside*
the engineered regions is a red flag: the graft perturbed the core.

**Radius of gyration** (`fig3_rg.png`) — flat means the trimer holds. Rising
means protomers separating. Within ~1–2% of the reference is fine.

**Per-region bars** (`fig4_region_bars.png`) — **the decision figure.** Mean
RMSF per region for each variant beside the reference. This is what belongs in
a thesis chapter.

**SASA** — a region that buries or exposes much more surface than the reference
has altered its antigenic presentation, which matters even if the fold is
stable.

**Inter-protomer H-bonds** — direct probe of interface integrity. A steady count
means the interface holds; a decline often precedes any movement in Rg.

### The decision rule

A candidate **passes** if all four hold:

1. RMSD plateaus by ~20 ns at a value comparable to the reference
2. Per-region RMSF not dramatically above the reference
3. Rg flat and within a couple of percent
4. Inter-protomer H-bond count stable

Failing any of these deprioritises the candidate **regardless of docking
score** — a variant that evades antibodies but falls apart is not useful.

---

## Limitations to state in your methods

- **Single 100 ns replicate** (unless you set `REPLICATES`) is a screen, not
  converged sampling. A difference comparable to run-to-run noise is not a real
  difference.
- **RMSF depends on the fitting reference and the discarded window.** Use the
  same `SKIP_NS` for every system and say what it was.
- **Every system must share the same preparation** — same force field, chain
  handling and protocol. LESSONS #7.
- **Rebuilt loops are models, not data.**
- **A walltime kill can leave a small gap between trajectory segments**, so
  frame spacing may not be perfectly uniform. This does not affect per-frame
  metrics.

---

## Resource sizing (ASPIRE2A, ~420k atoms, 1× A100)

| Phase | Wall time |
|---|---|
| EM | 10–20 min |
| NVT 100 ps | ~15 min GPU |
| Warm-up 50 ps | ~10 min GPU |
| NPT 1 ns | ~30–60 min GPU (**10 h on CPU** — this caused a walltime kill) |
| Production 100 ns | 40–60 GPU-h, ~63 ns/day, 2–3 × 24 h submissions |

Request `walltime=24:00:00`. 24 h is the maximum; 26 h is rejected.
