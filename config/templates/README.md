# Scaffold templates

A template is the full-length viral protein sequence with a single `{}`
placeholder marking the region to be humanised.

Example (HVR7 of HAdV-C5 hexon):

```
MATPSMMPQWSYMHISGQDASEYLSPGLVQFAR...{}...AEVKTHGKHYSYNSHW
```

Rules:
- Exactly one `{}` per template.
- Everything outside the placeholder is the conserved scaffold.
- The residue range the placeholder covers must match `region.range` in the
  config, otherwise downstream numbering (SASA, RMSF, AIRs) silently shifts.

Verify with:
```bash
python3 tool1_modelling/00_check_template.py config/hadv_c5_hvr7.yaml
```
