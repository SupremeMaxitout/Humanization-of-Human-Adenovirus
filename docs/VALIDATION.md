# Validation methodology

## What is checked, and why those things

Three categories, chosen because each predicts a specific downstream failure.

### Geometry — defects minimisation cannot repair

| Defect | Threshold | Why it is fatal |
|---|---|---|
| D-amino-acid centre | 0 | Inversion requires breaking bonds. Minimisation relaxes *around* it. Survives EM and restrained NVT, then explodes in free dynamics. |
| Heavy-atom overlap < 1.2 Å | 0 | The Lennard-Jones r⁻¹² term diverges; GROMACS reports `Epot = inf`. |
| Ring piercing | 0 | A bond threaded through an aromatic ring is topologically a knot. |
| Chain count ≠ 3 | fail | A homotrimer merged into one chain carries artificial peptide bonds across protomer junctions, stiffening the very interface an HVR study measures. |

Contacts of 1.8–2.3 Å between N/O atoms are **not** failures. They are ordinary
in predicted structures and resolve during energy minimisation. Treating them as
defects leads to "cleanup" steps that manufacture worse problems than they fix.

### Confidence — predictor self-assessment

| Metric | Default | Meaning |
|---|---|---|
| `mean_plddt` | ≥ 70 | Whole model. Below 70 is AlphaFold's own "low confidence" band. |
| `plddt_<region>` | ≥ 70 | **The decisive number.** A 90 global mean can hide a floppy graft. |
| `iptm` | ≥ 0.60 | Interface confidence — whether the trimer is predicted as an assembly. |
| `ptm` | ≥ 0.50 | Overall fold confidence. |
| `fraction_disordered` | ≤ 0.30 | AF3 estimate of disordered content. |
| `has_clash` | false | AF3's internal clash flag. |

pLDDT is read from the B-factor column, so AlphaFold3 and ColabFold are both
supported without special handling.

### SASA — is the graft still presented?

A humanised region only evades antibodies if it remains solvent-exposed.
Reported as a ratio to the wild-type reference. A ratio far below 1.0 means the
graft has folded inward and changed its antigenic presentation — structurally
fine, but potentially useless for the intended purpose.

This is comparative, not absolute. It informs interpretation rather than
gating, because there is no principled universal cutoff.

## Choosing thresholds

The defaults are conventional starting points, not derived constants. Two rules:

1. **Set them before looking at results.** Thresholds chosen after seeing the
   numbers turn a screen into post-hoc rationalisation. They live in the config
   and are written into the report's `Thresholds` sheet so the applied gate is
   recorded alongside the results.
2. **Include a negative control.** Run a deliberately bad graft — scrambled
   sequence, or a known-destabilising substitution — through the same gate. If
   it passes, the thresholds are too loose and the filter is not discriminating.

## Limitations

- Passing this gate means a model is **suitable for docking and MD**. It does
  not mean the design works.
- pLDDT is confidence, not accuracy. A confidently-predicted wrong structure
  scores well.
- SASA on a static model ignores loop dynamics; a region that looks buried may
  be transiently exposed. MD (Tool 4) is what addresses this.
- Region boundaries come from one numbering convention. A construct numbered
  differently silently misassigns every region — check `00_check_template.py`
  in Tool 1 first.
