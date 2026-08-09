# config

Run configuration. **This is the only file you edit to retarget the pipeline.**

| File | Purpose |
|---|---|
| `hadv_c5_hvr7.yaml` | the active run config |
| `templates/` | hexon scaffold sequences with one `{}` graft slot |

Every stage reads its region boundaries, paths and thresholds from here, so the
scripts cannot disagree about where a region starts or what counts as a pass.

## Key blocks

```yaml
target:
  virus: HAdV-C5
  protomer_length: 952        # folds chain B/C numbering onto protomer 1
  template: config/templates/hadv_c5_hexon_hvr7.txt

region:                       # the region being humanised
  name: HVR7
  query: <wild-type sequence>
  range: [418, 458]

regions:                      # ALL region boundaries - single source of truth
  HVR1: [136, 169]
  HVR7: [418, 458]

paths:
  models: data/models
  predictions: data/predictions

validation:
  reference: wild_type        # SASA baseline
  thresholds:                 # gate applied in stage 2
    min_mean_plddt: 70.0
    min_region_plddt: 70.0
    min_iptm: 0.60
```

Thresholds are set here, before results are seen, and copied into the output
report so the applied gate is recorded alongside the numbers.
