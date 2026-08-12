# Software versions

The exact versions used to produce the published results. Reproducing the
numbers requires these; other versions will usually work but may differ in
detail.

## Python environment

Pinned in [`requirements.txt`](requirements.txt) and
[`environment.yml`](environment.yml).

| Package | Version |
|---|---|
| Python | 3.10.9 |
| numpy | 1.26.4 |
| pandas | 2.2.2 |
| scipy | 1.13.1 |
| matplotlib | 3.8.4 |
| biopython | 1.83 |
| gemmi | 0.6.5 |
| openpyxl | 3.1.2 |
| PyYAML | 6.0.1 |

## External engines

| Software | Version | Notes |
|---|---|---|
| GROMACS | 2023.2 (GPU/CUDA 11.6) | module `gromacs/2023.2-gpu` on ASPIRE 2A |
| HADDOCK | 2.5, August 2025 release | personal licence from the Bonvin lab |
| CNS | 1.3, recompiled with HADDOCK 2.5 patches | free for non-profit use |
| AlphaFold3 | server (web), accessed 2026 | manual upload; no submission API |
| ColabFold | 1.5.5 (localcolabfold) | alternative, fully scriptable |

## Force field

**CHARMM27 + TIP3P** as shipped with GROMACS 2023.2.

Note: CHARMM36m is generally preferable for proteins with long disordered
loops — which is what an HVR is — but it is a separate port from the MacKerell
lab and must be installed into `$GMXLIB` manually. The results here used
CHARMM27, and `FORCE_FIELD=charmm36m` switches to CHARMM36m if you have
installed it. Whichever you use, apply it consistently to every system in a
comparison.

## Hardware

| | |
|---|---|
| Cluster | NSCC ASPIRE 2A |
| Scheduler | PBS Pro |
| GPU | NVIDIA A100-SXM4-40GB |
| CPU | AMD EPYC 7713 (16 cores per job) |
| Production throughput | ~63 ns/day per system (~420,000 atoms) |

## System size

| | |
|---|---|
| Protomer | 952 residues |
| Trimer | 3 chains, ~22,800 heavy atoms |
| Solvated | ~420,000 atoms (dodecahedral box, 1.5 nm padding) |
| Salt | 150 mM NaCl, neutralised |
| Temperature | 310 K |
