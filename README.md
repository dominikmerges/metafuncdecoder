# MetaFuncDecoder

MetaFuncDecoder converts raw functional annotations (KO terms, Pfam IDs, CAZy families)
from JGI KOfamScan, HMMER Pfam, and dbCAN files into ecologically interpretable
functional profiles with cross-database confidence scoring.

## What it does

Many (meta)genome annotation pipelines produce files full of identifiers
(K03578, PF00001, GH18) that are hard to interpret directly. MetaFuncDecoder:

1. Loads annotation files from JGI, HMMER, dbCAN, or a pre-merged table
2. Maps identifiers to functional categories defined in user-editable YAML files
3. Scores confidence based on how many databases agree on the same category
4. Outputs a single CSV with one row per gene per matched category, including subcategory labels and evidence from each database

Functional categories are defined in plain YAML files that can be copied and
modified without touching any Python code. Default categories (carbon cycling,
nitrogen cycling, antibiotic resistance) include citations for their pattern
definitions.

## Installation

Clone or download this repository, then install dependencies:

```bash
pip install -r requirements.txt
```

For development and running the test suite, also install:

```bash
pip install -r requirements-dev.txt
```

Then download the reference databases:

```bash
bash download_dbs.sh
```

This downloads three files into `data/`:
- `pfam2go.txt` — Pfam to GO term mapping (Gene Ontology Consortium)
- `FamInfo.txt.04232020.tsv` — CAZy family activity descriptions (dbCAN, April 2020)
- `ko_list` — KEGG ortholog definitions (KOfam; free for academic use)

## Quick start

**Join mode** — separate annotation files, joined on gene ID:

```bash
# KO file only — single category
python metafuncdecoder.py \
  --ko-file 80425.assembled.ko \
  --function carbon_cycling

# KO + Pfam — all categories in one run
python metafuncdecoder.py \
  --ko-file 80425.assembled.ko \
  --pfam-file 80425.assembled.pfam.blout \
  --function all \
  --output results/

# KO + Pfam + CAZy — select specific categories
python metafuncdecoder.py \
  --ko-file 80425.assembled.ko \
  --pfam-file 80425.assembled.pfam.blout \
  --cazy-file dbcan_overview.txt \
  --function carbon_cycling nitrogen_cycling \
  --output results/

# Conservative run — require 2 dbCAN tools for CAZy; exclude broad-evidence hits
python metafuncdecoder.py \
  --ko-file 80425.assembled.ko \
  --pfam-file 80425.assembled.pfam.blout \
  --cazy-file dbcan_overview.txt \
  --function all \
  --cazy-min-tools 2 \
  --exclude-broad \
  --output results/
```

**Combined mode** — single pre-merged annotation table:

```bash
python metafuncdecoder.py \
  --combined-table annotations.tsv \
  --ko-col ko_id \
  --pfam-col pfam_hits \
  --cazy-col cazy_id \
  --function all
```

## Input file formats

### KO file (`--ko-file`)
JGI KofamScan output (`.assembled.KO`). Tab-separated, no header.
Only rows with `Yes` in the threshold column are used; `No` rows are discarded.
KofamScan applies a per-model score cutoff derived from the KOfam HMM profiles,
so the `Yes`/`No` column already reflects a model-specific significance decision.
No additional score filtering is needed or applied.

```
Ga0105250_100000012  Yes  KO:K04775  94.31  1  281  1  281  6.3e-200  624  281
```

### Pfam file (`--pfam-file`)
Two formats are supported, auto-detected:
- **JGI blout** (`.assembled.pfam.blout`): tab-separated, gene_id first, pfam_id second.
  Pfam IDs in `pfamXXXXX` format are normalised to `PFXXXXX` automatically.
- **HMMER domtblout**: standard `--domtblout` output; comment lines start with `#`.

MetaFuncDecoder uses all Pfam hits present in the file without applying a score
filter. Whether hits are pre-filtered depends on how HMMER was run upstream:

- **Recommended:** run HMMER with `--cut_ga` (Pfam gathering cutoffs). This applies
  a per-model significance threshold equivalent to what KOfamScan does for KEGG, and
  all hits in the output file can be treated as reliable.
- **Without `--cut_ga`:** the file may contain low-scoring hits. Pre-filter the blout
  or domtblout file before passing it to MetaFuncDecoder, or treat Pfam-only
  (`confidence_level = low`) calls with additional caution.

JGI annotation pipelines typically apply gathering cutoffs, so JGI blout files are
generally pre-filtered. If you are unsure, check the HMMER command used to generate
your Pfam file.

### CAZy file (`--cazy-file`)
dbCAN `overview.txt` output. Tab-separated with header.
The CAZy family column is used for family assignments. MetaFuncDecoder accepts
both column names produced by different dbCAN versions: `HMMER` and `dbCAN`.

```
Gene_ID   EC#   HMMER          dbCAN_sub   DIAMOND   #ofTools
gene_001  3.2   GH18(19-390)   GH18_e1     GH18      3
```

### Combined table (`--combined-table`)
Any tab-separated table with annotation columns. Use `--ko-col`, `--pfam-col`,
`--cazy-col`, and `--gene-id-col` to specify the relevant column names.

By default, MetaFuncDecoder exits with an error if a specified annotation column is
not found in the table. This prevents silent runs where a missing column is skipped without warning. To allow missing columns (and skip the corresponding
database with a warning), pass `--allow-missing-cols`:

```bash
python metafuncdecoder.py \
  --combined-table annotations.tsv \
  --ko-col ko_id \
  --pfam-col pfam_hits \
  --cazy-col cazy_id \
  --function all \
  --allow-missing-cols
```

## Output files

### Console output

MetaFuncDecoder always prints a **Functional Annotation Summary Report** to stdout:

```
======================================================================
FUNCTIONAL ANNOTATION SUMMARY REPORT
======================================================================

Function: carbon_cycling
Total annotations: 5683
  metagenome: 100.0% (5683 annotations)
```

### CSV output

All CSV outputs go to `--output` (default: `results/`).
The output file is named `<dataset>_confidence_annotations.csv`, where `<dataset>` is
`metagenome` by default and can be set with `--dataset-name`. One row is written per gene per
functional category matched; genes with no match are omitted by default (use
`--include-unannotated` to retain them).

### Output columns

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier (from input files) |
| `functional_category` | Broad functional category matched (e.g. `carbon_cycling`) |
| `subcategories` | Functional subcategories, semicolon-separated. For CAZy hits: enzyme class (e.g. `chitinase;cellulase`). For KEGG/Pfam hits: process group if the matched pattern carries a subcategory label (e.g. `nitrogen_fixation`, `antibiotic_inactivation`). **May be empty** for genes matched only by broad patterns (e.g. `GO:0005975` carbohydrate metabolic process, or non-specific text patterns) that do not resolve to a specific enzyme class. See `evidence_specificity`. |
| `evidence_specificity` | `specific` if at least one database matched a substrate-specific pattern (named enzyme KO term, CAZy family with defined activity, or Pfam GO term naming the function directly) and a subcategory label was assigned. `broad` if the category match came only from non-specific patterns (e.g. `GO:0005975`, `antibiotic`, `antimicrobial`) with no subcategory. Empty string for unannotated genes. Use `--exclude-broad` to omit broad hits from the output entirely, or filter post-hoc on `evidence_specificity != 'broad'`. |
| `subcategory_n_databases` | Maximum number of databases that independently support the same subcategory for this gene (0–3). A value of 2 means at least two databases (e.g. KEGG and CAZy) assigned the same subcategory label (e.g. `chitinase`). Complements `n_databases`, which reflects broad-category agreement: `n_databases = 2, subcategory_n_databases = 0` means two databases agree on the broad category but support different or no subcategories. |
| `cazy_families` | CAZy family IDs from the dbCAN HMMER column (e.g. `GH18;GT2`) |
| `cazy_subfamilies` | CAZy subfamily IDs from the dbCAN `dbCAN_sub` column (e.g. `GH5_21`), semicolon-separated. Only populated when `--cazy-file` contains a `dbCAN_sub` column. When present, subfamily-level lookup is used for functional assignment; family-level is used as fallback. |
| `cazy_ec` | EC numbers from the dbCAN `EC#` column, semicolon-separated. Only populated when using `--cazy-file`. |
| `cazy_tool_count` | Number of dbCAN tools agreeing (1–3). Only populated when using `--cazy-file`. dbCAN recommends ≥ 2 for reliable assignments. |
| `pfam_ids` | Pfam accession numbers (e.g. `PF00128`) |
| `ko_terms` | KEGG Orthology terms (e.g. `K03791`) |
| `cazy_description` | Activity description(s) from CAZy FamInfo, pipe-separated |
| `pfam_go_names` | GO term name(s) mapped via Pfam2GO, pipe-separated |
| `ko_definition` | KO definition string from the KEGG `ko_list` |
| `supporting_databases` | Which databases provided evidence, semicolon-separated (`CAZy`, `Pfam`, `KEGG`) |
| `n_databases` | Number of databases that independently support the annotation (1–3) |
| `confidence_level` | `high` / `medium` / `low` / `none` (see table below) |
| `confidence_score` | Numeric confidence: 1.0 / 0.67 / 0.33 / 0.0 |

### Confidence scoring

| Level | Meaning | Score |
|-------|---------|-------|
| high | All 3 databases agree | 1.0 |
| medium | 2 databases agree | 0.67 |
| low | 1 database | 0.33 |
| none | No functional annotation | 0.0 |

### Reproducibility and provenance

Every run writes a `run_manifest.json` alongside the output CSV. This file records:
- MetaFuncDecoder version
- Full CLI arguments
- Selected functional categories
- Paths and SHA-256 checksums of all reference database files and YAML category files

Keep `run_manifest.json` alongside the result CSV when archiving data for publication.
It allows any run to be exactly reproduced or audited without relying on memory of
which database version was used.

### Interpreting low-confidence hits

`low` confidence (single-database) hits include two distinct types of evidence,
distinguished by the `evidence_specificity` column:

- **`specific`**: the gene matched a substrate-specific pattern — a named enzyme KO
  term, a CAZy family with a defined activity, or a Pfam GO term that directly names
  the function. A subcategory label is present. These are reliable single-database calls.
- **`broad`**: the gene matched only a non-specific pattern — a broad GO term
  (`GO:0005975` carbohydrate metabolic process), or a text pattern such as `antibiotic`
  or `antimicrobial`. No subcategory label is assigned. These hits indicate a gene is
  plausibly associated with the broad functional category but the evidence does not
  support a specific enzyme class. Treat with additional caution in downstream analysis.

To exclude broad hits:

```bash
# At run time — broad rows never written to the output CSV
python metafuncdecoder.py --ko-file my.ko --function all --exclude-broad

# Post-hoc in Python
df = df[df["evidence_specificity"] != "broad"]

# Or keep only medium/high confidence (which also excludes all broad hits,
# since broad hits are always single-database and therefore always low confidence)
df = df[df["confidence_level"].isin(["medium", "high"])]
```

## Functional categories

Categories are defined in YAML files in `categories/`. Three defaults are provided:

| File | Category | Basis |
|------|----------|-------|
| `carbon_cycling.yaml` | Carbohydrate-active enzyme functions | López-Mondéjar et al. (2022) mSystems; CAZy database |
| `nitrogen_cycling.yaml` | Nitrogen transformation functions | NCycDB (Tu et al. 2019) gene families |
| `antibiotics.yaml` | Antibiotic resistance and biosynthesis | CARD/ARO (Alcock et al. 2023) |

Each YAML file contains separate pattern lists for KEGG definitions, Pfam GO
term names, and CAZy activity descriptions, plus exact GO term IDs for Pfam
matching. See `categories/README.md` for the full format.

### Running multiple categories

`--function` accepts one or more category names, or `all` to run every YAML in
`categories/` in a single pass:

```bash
# All categories — one CSV with a functional_category column to distinguish them
python metafuncdecoder.py --ko-file my.ko --function all

# Selected categories
python metafuncdecoder.py --ko-file my.ko --function carbon_cycling nitrogen_cycling
```

The output CSV contains only the rows for the requested categories. Use
`--function all` to include every loaded category.

### Adding a custom category

```bash
cp categories/carbon_cycling.yaml categories/my_function.yaml
# edit name, description, citations, and patterns
python metafuncdecoder.py --ko-file my.ko --function my_function
```

### Custom pattern search (without editing YAML)

```bash
python metafuncdecoder.py \
  --ko-file my.ko \
  --custom-ko-pattern "laccase" \
  --custom-category-name "laccase_activity"
```

## Reference databases

```bash
bash download_dbs.sh                        # download to data/ (default)
bash download_dbs.sh --data-dir /my/path/   # custom location
```

| Database | Source | License |
|----------|--------|---------|
| Pfam2GO | Gene Ontology Consortium | CC BY 4.0 |
| CAZy FamInfo | dbCAN (Yin et al. 2012) | Open access |
| KEGG ko_list | KOfam (Aramaki et al. 2020) | Free for academic use |

## Development

To run the test suite:

```bash
pytest -q
```

All 66 tests run in-memory and do not require the reference databases to be downloaded.

## Citation

If you use MetaFuncDecoder, please cite the software release:

> Merges D (2026). MetaFuncDecoder v1.0.0. Zenodo. https://doi.org/10.5281/zenodo.19009291

And the relevant category definition sources:
- López-Mondéjar et al. (2022) mSystems 7:e00829-22
- Tu et al. (2019) Bioinformatics 35:1040-1048
- Alcock et al. (2023) Nucleic Acids Res. 51:D690-D699
- Aramaki et al. (2020) Bioinformatics 36:2251-2252 (KOfam)
- Lombard et al. (2014) Nucleic Acids Res. 42:D490-D495 (CAZy)
