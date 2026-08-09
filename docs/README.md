# docs

Methodology and hard-won operational knowledge.

| Doc | Covers |
|---|---|
| `VALIDATION.md` | what stage 2 checks, why those things, how to choose thresholds |
| `LESSONS.md` | every failure encountered building this, and its fix |
| `TROUBLESHOOTING.md` | symptom to cause to fix; when to resume vs restart |
| `PARAMETERS.md` | which MD parameters to change, which never to touch |
| `STRUCTURE_PREP.md` | cleaning structures, PyMOL recipes, missing residues |
| `ANALYSIS.md` | what each MD output means and how to read it |

`LESSONS.md` is the one to read first if you are reproducing this work. It
documents roughly two weeks and 15,000 CPU-hours of failures - a missing
`refcoord-scaling = com`, a conjugate-gradient stage that could not survive
water constraints, and a "cleanup" step that manufactured the defects it was
meant to remove. Every one is now either fixed in code or checked automatically.
