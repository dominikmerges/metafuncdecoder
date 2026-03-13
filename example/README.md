# Example input files

Small subset of a real JGI metagenome annotation (IMG taxon 80425, soil metagenome)
for quick-start testing. 520 genes, two annotation files.

## Files

| File | Format | Description |
|------|--------|-------------|
| `example.ko` | JGI KofamScan `.assembled.KO` | KEGG ortholog assignments |
| `example.pfam.blout` | JGI HMMER blout | Pfam domain hits |

The subset includes genes annotated to carbon cycling (cellulase, chitinase, xylanase),
nitrogen cycling (nitrogenase, nitrate reductase), and antibiotic resistance functions,
alongside unannotated genes — representative of a real metagenome.

## Quick start

From the repository root, after running `bash download_dbs.sh`:

```bash
# KO only
python metafuncdecoder.py \
  --ko-file example/example.ko \
  --function carbon_cycling

# KO + Pfam — all categories
python metafuncdecoder.py \
  --ko-file example/example.ko \
  --pfam-file example/example.pfam.blout \
  --function all \
  --output results/example/

# Single category
python metafuncdecoder.py \
  --ko-file example/example.ko \
  --function carbon_cycling
```

Expected output: `results/example/metagenome_confidence_annotations.csv`

## Source

Derived from JGI IMG metagenome 80425 (public dataset). Full annotation files
available via the JGI IMG/M portal.
