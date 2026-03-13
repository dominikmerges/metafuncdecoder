# Changelog

All notable changes to MetaFuncDecoder are documented here.

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
