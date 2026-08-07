# data/

Working directory. **Nothing here is committed** (see `.gitignore`) — the
proteome FASTA is ~110 MB uncompressed and predictions run to gigabytes.

## Getting the human proteome

```bash
cd data
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
```

Then `01_parse_proteome.py` converts it to `ensembl_proteins.json`.

Record the Ensembl release you used — it goes in the manifest and belongs in
your methods, because candidate hits depend on the proteome version.

## Layout created by the pipeline

```
data/
  Homo_sapiens.GRCh38.pep.all.fa   you download this
  ensembl_proteins.json            step 01
  predictions_raw/                 AF3 zips, or colabfold_batch output
  predictions_extracted/           step 04
  best_models/                     step 05  (mmCIF)
  trimers/                         step 06  (PDB, 3 chains)  -> MD
  trimers_docking/                 step 07  (PDB, 1 chain)   -> HADDOCK
```
