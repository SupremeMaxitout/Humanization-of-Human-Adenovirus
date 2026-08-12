# Lessons learned

The MD stage is the second version of this pipeline. The first took roughly two
weeks and ~15,000 CPU-hours to get three of six systems through equilibration,
and almost none of that time was spent on science. Every failure below is now
either fixed in code or checked automatically by the preflight in `setup_md.sh`.

If you are reproducing this work, reading this file will save you those weeks.

---

## #1 — THE decisive bug: `refcoord-scaling` (cost ~2 weeks)

**Symptom.** Energy minimisation fine. NVT fine (100 ps, no complaints). Then
NPT dies between step ~800 and ~1500:

```
Step 1140  LINCS WARNING
relative constraint deviation after LINCS:
rms 14898.8, max 2061565.1 (between atoms 7261 and 7262)
bonds that rotated more than 30 degrees
   7261 7262  90.0  0.1773  222649.1406  0.1080
Segmentation fault (core dumped)
```

One bond stretching to 222,649 nm while the rest of the system is perfectly
stable (rms deviation 0.00003).

**Cause.** `npt.mdp` had `define = -DPOSRES` and `pcoupl = C-rescale` but **no
`refcoord-scaling` line**. GROMACS then defaults to `No`. The barostat rescales
every atom's coordinates as the box changes, while the position-restraint
*reference* coordinates stay frozen at their original values. Each step the
restrained atoms are pulled further toward references that no longer correspond
to anything. Tension accumulates for ~1 ps, then one bond gives way.

**Fix.** One line in every restrained stage that also has a barostat:

```
refcoord-scaling = com
```

**Why it took two weeks.** The failure *looks* exactly like a bad structure: a
single localised explosion. We chased D-amino acids, ring piercings, steric
clashes, barostat choice, timestep and force field. The tell we missed: NVT (no
barostat) always passed, and a 50 ps warm-up that happened to include
`refcoord-scaling = com` also passed. The same geometry surviving one stage and
dying in the next should have pointed straight at the mdp difference.

**Now enforced.** `setup_md.sh` refuses to start if any restrained NPT stage has
a barostat without `refcoord-scaling = com`.

> A note on the folk explanation. It is tempting to summarise this as "we just
> used the AlphaFold output directly and it worked". That is half true. Using
> raw predicted output removed damage we had introduced ourselves (see #3), but
> systems built from raw output *still* died at ~1400 steps until
> `refcoord-scaling = com` was added. Both changes were needed; the mdp line
> was decisive. Record it that way in your methods.

---

## #2 — Never add a conjugate-gradient minimisation stage

A second-stage `integrator = cg` minimisation was added to "clean up residual
clashes". It produced:

```
Fatal error:
The coordinates could not be constrained.
Minimizer 'cg' can not handle constraint failures,
use minimizer 'steep' before using 'cg'.
```

TIP3P water is held rigid by SETTLE regardless of `constraints = none`, and
conjugate gradient cannot recover when a water constraint fails during its line
search. In a solvated system, `cg` is inherently fragile.

Worse, the stage was **also useless**: NVT read `em.gro` (the steepest-descent
output), not the cg output, so its result was discarded even when it succeeded.

---

## #3 — Do NOT "pre-relax" model structures in vacuum

Input PDBs were run through an OpenMM minimisation with
`nonbondedMethod=NoCutoff` and no solvent, intending to clean up clashes. In
vacuum there is no water to screen electrostatics, so every charged/polar pair
(Arg/Lys to backbone carbonyl, Asp, Glu) attracts with full unscreened force.
The minimiser dragged them into contacts as short as **1.77 Å** where a real
salt bridge sits near 2.8 Å.

Measured on the same file:

| | raw AlphaFold3 | after vacuum relax |
|---|---|---|
| D-amino-acid centres | 0 | 8 |
| hard clashes < 2.3 Å | 6 (mild, 1.8–2.0 Å) | 31 (down to 1.77 Å) |

**The preparation created the damage it was meant to remove.** Those
over-compressed pairs are held in place by restraints through EM and NVT, then
recoil violently once NPT gives them freedom.

**Fix.** Feed model coordinates straight into `pdb2gmx`. GROMACS's own EM, in
explicit solvent with PME, resolves mild model contacts perfectly well.

Related trap: `constraints=HBonds` in an OpenMM minimisation of a strained
structure makes the energy *increase* — constraint enforcement fights the
minimisation. Two of three systems came out at 2.7e7 kJ/mol, unchanged from
their start, while appearing to have "run successfully".

---

## #4 — A homotrimer is three chains, not one

The first pipeline used `gmx pdb2gmx -merge all`, which fuses all three
protomers into a single molecule with **artificial peptide bonds between
chains**. Those fake covalent links stiffen precisely the inter-protomer
interface that an HVR flexibility study is measuring.

**Fix.** `-merge no`. You should see `topol_Protein_chain_A.itp`, `_B`, `_C`
after `pdb2gmx`. If you see only one, stop.

---

## #5 — PBS and GROMACS must not share a log filename

`job_prod.pbs` wrote PBS stdout to `prod.log` — the same file GROMACS uses. PBS
overwrote it, and on the next resume GROMACS refused to continue:

```
Can't read 1048576 bytes of 'prod.log' to compute checksum.
The file has been replaced or its contents have been modified.
Cannot do appending because of this condition.
```

The job exits in seconds. It *looks* like it ran; nothing advances. Two full
queue cycles were lost to this before it was spotted.

**Fix.** PBS writes to `prod_pbs.log`, and production runs with `-noappend` so
each window writes `prod.partNNNN.*` instead of appending.

---

## #6 — Know what "finished" actually means

`Finished mdrun on rank 0` prints after a `-maxh` graceful stop as well as after
a genuine completion. It is **not** a completion marker, and treating it as one
led to a system being called done at 97.5 ns.

**The real test is both of:**
- `prod*.gro` exists (the final `.gro` is only written when `nsteps` is reached)
- the last logged step equals `nsteps`

`tools/check_status.sh` applies both.

---

## #7 — Prepare every system in the comparison identically

Three systems were built with `-merge all` and an NPT lacking
`refcoord-scaling`; three were rebuilt properly. Comparing RMSF across those two
preparations is a methods flaw a reviewer will catch, because the old three have
artificial inter-chain bonds exactly where flexibility is being measured.

**Rule.** If you change the preparation for one system, rebuild all of them.
Cheap now that setup runs on GPU in 1–2 h; expensive after you have published.

---

## #8 — Distinguish "crashed" from "not finished yet"

Time was lost to a status check that only tested for `npt.gro`:

- Jobs still sitting in the queue were reported as `[FAIL]`.
- `stepNNNb.pdb` crash dumps from an *old* attempt survived in the run
  directory and made a healthy run look like a fresh explosion.
- A job killed at step 480,000 of 500,000 (96% done, checkpoint intact) was read
  as a failure rather than "resubmit and it finishes in 20 minutes".

**Fix.** `tools/check_status.sh` reports the actual stage, queue state and ns
completed. `tools/diagnose.sh` separates a walltime kill from a real crash
before you touch anything.

---

## #9 — Analysis without PBC correction is wrong

The original `analyse_md.sh` ran `gmx rms` directly on `prod.xtc`. Over 100 ns a
trimer diffuses across the periodic boundary and protomers get wrapped to
opposite sides of the box, producing large fake jumps in RMSD and Rg.

**Fix.** `gmx trjconv -pbc whole -center` then `-fit rot+trans` before computing
any metric. Built into the current `analyse_md.sh`.

And with `-noappend`, analysis must use the **joined** trajectory
(`tools/concat_trajectories.sh`). Otherwise you silently analyse 60 ns of one
system against 100 ns of another.

---

## #10 — Cluster habits that cost real time

| Mistake | Consequence | Fix |
|---|---|---|
| Running compute on the login node | Formal warning; processes killed; 3 strikes can block the account | `qsub -I` for interactive, `qsub` for batch. Preflight now refuses. |
| `qsub -q g1` | `Access to queue is denied` | Submit to `normal`; PBS routes GPU jobs into `g1` |
| Requesting 26 h walltime | Rejected — 24 h is the maximum | `-maxh 23.0` + `-cpi`; 2–3 submissions is the designed path |
| 10 h walltime with CPU-only equilibration | NPT killed at 96% | 24 h + GPU offload → setup takes 1–2 h |
| `rm -rf` the directory you are standing in | `getcwd: No such file or directory`; `qsub` stops working | `cd ~` first |
| Long job in an SSH session | VPN drop kills it | Batch jobs, or `nohup … &` |
| `scp`-ing a file and not verifying | A missing `npt_warmup.mdp` wasted a full queue cycle | Preflight checks all five mdp files |
| Filename with `.cif` in the middle | PDBFixer parses `x_model_0.cif_single_chain.pdb` as mmCIF and crashes | Short, clean filenames |
| `from openmm.unit import *` | Shadows Python's `sum`; `TypeError: object of type 'generator' has no len()` | Import explicitly |

---

## The five checks that would have prevented all of it

1. `grep refcoord-scaling npt.mdp` whenever `-DPOSRES` meets a barostat.
2. Run a geometry check on the input **before** queuing anything.
3. Never modify coordinates in vacuum.
4. Confirm `qstat` is empty before calling anything a failure.
5. Check file timestamps before trusting `step*.pdb` crash dumps.
