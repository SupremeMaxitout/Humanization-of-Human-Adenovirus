# Depositing: DOIs for the code and the data

Two separate deposits, for two different things:

1. **The code** — this repository, archived from GitHub. Gives you a citable
   software DOI.
2. **The data** — trajectories, models and docking output. Far too large for
   git, but without them nobody can verify your analysis.

Both are free. Neither needs institutional permission.

---

## 1. Code DOI (Zenodo ↔ GitHub)

Takes about ten minutes, once.

**Link the accounts**
1. Sign in at [zenodo.org](https://zenodo.org) with your GitHub account.
2. Go to **Settings → GitHub**.
3. Find `Humanization-of-Human-Adenovirus` and toggle it **ON**.

Zenodo now watches the repository. Nothing is archived until you cut a release.

**Cut the release**
1. On GitHub: **Releases → Create a new release**.
2. Tag `v1.0.0`, title `v1.0.0 — complete four-stage pipeline`.
3. In the description, summarise what the release contains and what it was used
   for. Mention the four stages and the system it was developed on.
4. **Publish release.**

Within a minute or two Zenodo archives that snapshot and mints a DOI.

**Add the badge.** Zenodo shows a DOI badge on the deposit page. Put its
markdown at the top of the README, and record the DOI in `CITATION.cff`.

**Versioning.** Every subsequent release gets its own DOI, plus a *concept DOI*
that always resolves to the latest. Cite the concept DOI in a CV; cite the
specific version DOI in a paper.

---

## 2. Data DOI (Zenodo dataset)

Trajectories are the actual product of the MD stage and cannot go in git. A
100 ns trajectory of a ~420,000-atom system is roughly 0.5–1 GB compressed.

Zenodo accepts 50 GB per record by default, which is enough for six systems if
you deposit the compressed trajectories and not the raw checkpoints.

**What to deposit**

| Include | Why |
|---|---|
| `prod_full.xtc` per system | the trajectory itself |
| `prod.tpr` per system | required to read the trajectory |
| final `.gro` per system | last frame |
| all `.mdp` files | the exact parameters used |
| `topol.top` + `*.itp` | the topology |
| `analysis/*.xvg` and figures | derived results |
| input PDBs | starting structures |

| Exclude | Why |
|---|---|
| `.cpt`, `.trr` | huge, only useful for restarting |
| `work_pbs/` | HADDOCK scratch |
| `.log` files | verbose, low value |

**Package it**

```bash
cd ~/md_runs
for s in wild_type fold_013 fold_022 fold_039 fold_042 fold_044; do
  mkdir -p ~/deposit/$s
  cp $s/prod_full.xtc $s/prod.tpr $s/topol.top ~/deposit/$s/
  cp $s/*.itp $s/*.mdp ~/deposit/$s/ 2>/dev/null
  ls $s/prod*.gro >/dev/null 2>&1 && cp $s/prod*.gro ~/deposit/$s/
done
cp -r ~/md_runs/analysis ~/deposit/
cd ~/deposit && tar czf hadv_md_trajectories_v1.tar.gz */ && ls -lh *.tar.gz
```

Include a `README.txt` in the archive describing each directory, the force
field, the simulation length, and the software versions.

**Upload**
1. [zenodo.org/uploads/new](https://zenodo.org/uploads/new)
2. Upload type: **Dataset**
3. Title: `Molecular dynamics trajectories: HAdV-C5 hexon HVR humanisation variants`
4. Description: systems, length, force field, conditions, software versions
5. Under **Related identifiers**, add the code DOI as *"is supplemented by"*
6. Licence: CC-BY-4.0 is standard for research data
7. Publish

Then link the data DOI from the README so the two records point at each other.

---

## 3. Citing it

**On a CV, under Software / Research Outputs:**

> Khan, M. (2026). *Humanisation of Human Adenovirus: a reproducible pipeline
> for immune-evasive capsid engineering* (v1.0.0) [Computer software].
> Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

**In a thesis or paper methods section:**

> Structure preparation, validation, docking and molecular dynamics were
> performed using the Humanisation of Human Adenovirus pipeline v1.0.0
> (DOI: 10.5281/zenodo.XXXXXXX). Trajectories are deposited at
> DOI: 10.5281/zenodo.YYYYYYY.

A DOI is what turns a GitHub repository into a citable research output. It is
permanent, indexed, and does not disappear if the repository is renamed or
deleted — which is exactly why reviewers and hiring panels take it more
seriously than a bare link.

---

## 4. Optional: ORCID

If you do not have an [ORCID](https://orcid.org) yet, get one — it is free and
takes two minutes. Add it to your Zenodo profile and to `CITATION.cff`, and
every deposit is automatically attached to your researcher identity.
