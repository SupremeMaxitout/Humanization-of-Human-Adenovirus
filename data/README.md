# data

Working directory. **Contents are git-ignored** - nothing here is version
controlled except this file.

Large inputs and intermediates live here so the repository stays small and
clonable.

## Expected layout

```
data/
  Homo_sapiens.GRCh38.pep.all.fa   human proteome (download, ~110 MB)
  ensembl_proteins.json            parsed proteome (stage 1)
  predictions/                     raw AlphaFold3 / ColabFold output
  models/                          selected models, PDB, multi-chain
  models_docking/                  same models merged to single chain
```

## Getting the proteome

```bash
cd data
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
```

## Why nothing is committed here

Trajectories, proteomes and prediction archives are far too large for git and
would make the repository unusable. For published work, deposit these on Zenodo
and cite the DOI.
