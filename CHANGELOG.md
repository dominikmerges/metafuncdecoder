# Changelog

All notable changes to MetaFuncDecoder are documented here.

## [2.0.0] — 2026-10-27

### Added

- **`--viral`**: adds a `viral_cazy_role` column labelling each CAZy family by its known phage role (`peptidoglycanase_lysis`, `virion_structural_or_lysis`, `polysaccharide_degrading`, `non_catalytic_cbm`, `other_cazy`). Family sets fixed to Merges et al. 2026 AEM, Table S15. The role is per gene and uses all of the gene's CAZy hits. Off by default; output is unchanged without the flag.

### Removed

- **`confidence_level` and `confidence_score` columns.** Use `n_databases` and `subcategory_n_databases` instead. Breaking change for scripts that read these columns.

### Changed

- **Output file renamed** from `<dataset>_confidence_annotations.csv` to `<dataset>_annotations.csv`. Breaking change for pipelines that expect the old name.

### Fixed

- `carbon_cycling.yaml`: CE9 was inside the GH153 comment and so missing from the `peptidoglycanase` subcategory; it is now listed, matching Merges et al. 2026 AEM, Table S15.
- Combined-table mode now splits `+`-joined CAZy domains (e.g. `CBM50+GH24`) and strips coordinates, as `--cazy-file` already did.
- dbCAN subfamilies are matched to their family by prefix. Before, they were matched by position, so a gene like `CBM2;GH5` with subfamily `GH5_7` could attach `GH5_7` to `CBM2`.
- With `--include-unannotated`, unannotated rows now get `subcategory_n_databases = 0`. Before, the column became decimals with blanks (e.g. `1.0`) and moved to the end of the table when the first gene was unannotated.

## [1.0.0] — 2026-03-13

Initial public release.

### Added

- **Join mode**: separate `--ko-file`, `--pfam-file`, and `--cazy-file` inputs, auto-joined on gene ID
- **Combined mode**: pre-merged annotation table via `--combined-table` with configurable column names
- **YAML functional categories**: three default categories (`carbon_cycling`, `nitrogen_cycling`, `antibiotics`); user-editable without modifying Python code
- **Three-database confidence scoring**: `high` / `medium` / `low` / `none` based on how many databases (KEGG, Pfam, CAZy) agree on the same functional category
- **`evidence_specificity` column**: distinguishes substrate-specific matches (`specific`) from broad/non-specific matches (`broad`); empty for unannotated genes
- **`subcategory_n_databases` column**: maximum number of databases independently supporting the same subcategory label (complements `n_databases`)
- **`--cazy-min-tools N`**: per-family filter requiring at least N dbCAN tools to agree before a CAZy hit contributes evidence; handles both integer (`"3"`) and float-formatted (`"3.0"`) tool counts from combined-mode tables
- **`--exclude-broad`**: removes broad-evidence rows from the output CSV at run time
- **`--include-unannotated`**: retains genes with no functional match in the output (excluded by default)
- **`--dataset-name`**: sets the output CSV filename prefix (default: `metagenome`)
- **`--allow-missing-cols`**: in combined mode, skip missing annotation columns with a warning instead of exiting with an error
- **`run_manifest.json`**: written alongside each output CSV; records script version, CLI arguments, selected categories, and SHA-256 checksums of all reference DB and YAML files
- **`download_dbs.sh`**: downloads all three reference databases into `data/`; auto-converts CAZy FamInfo XLS to TSV
- **CAZy subfamily support**: `dbCAN_sub` column parsed when present; subfamily-level lookup used with family-level fallback
- **Dual dbCAN column name support**: accepts CAZy family column named `HMMER` or `dbCAN` (both produced by different dbCAN versions)
- **`--custom-ko-pattern` / `--custom-category-name`**: ad hoc pattern search without editing YAML files
- 66 unit tests covering all core logic paths

### Known limitations

- Pfam score filtering is not applied internally; reliability of Pfam hits depends on whether HMMER was run with `--cut_ga` gathering cutoffs upstream
- `print` / `sys.exit` used for user-facing messages and fatal errors (library-style refactor deferred to v2)
- YAML schema validation is pattern-level only; file-level hard-fail for malformed YAML deferred to v2
