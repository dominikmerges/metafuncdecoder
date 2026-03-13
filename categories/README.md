# Functional Category Files

Each `.yaml` file in this directory defines one functional category used to
interpret annotations from KEGG, Pfam (via GO terms), and CAZy.

## Structure

```yaml
name: category_name
description: "Human-readable description"
citations:
  - "Author et al. (year) Journal"

kegg_patterns:       # regex matched against KEGG definition strings
  - "simple pattern"                     # plain string — no subcategory label
  - pattern: "grouped pattern"           # dict form — adds a subcategory label
    subcategory: subcategory_name        # written to the 'subcategories' output column

pfam_patterns:       # regex matched against GO term name strings (Pfam2GO)
  - "simple pattern"
  - pattern: "grouped pattern"
    subcategory: subcategory_name

pfam_go_terms:       # exact GO ID matches (Pfam2GO); no subcategory labels here
  - "GO:XXXXXXX"    # comment describing the term

cazy_subcategories:  # CAZy-specific subcategories (roll up to this broad category)
  subcategory_name:
    description: "..."
    patterns:        # matched against CAZy FamInfo activity description strings
      - "pattern"
```

All regex patterns are case-insensitive and follow Python `re` syntax.

### Subcategory labels

Both `kegg_patterns` and `pfam_patterns` support an optional `subcategory` field.
When a pattern with a subcategory label matches, that label is written to the
`subcategories` column of the output alongside the broad `functional_category`.

This allows finer-grained grouping without creating separate YAML files. For
example, nitrogen cycling patterns are grouped into `nitrogen_fixation`,
`nitrification`, `denitrification`, `dnra`, `anammox`, and `ammonia_assimilation`.
Antibiotic patterns are grouped into `antibiotic_inactivation`,
`antibiotic_target_alteration`, `antibiotic_efflux`, and `antibiotic_biosynthesis`.

Plain string patterns (no subcategory) still match the broad category but leave
`subcategories` empty for that hit.

## Provided categories

| File | Category | Databases used |
|------|----------|----------------|
| `carbon_cycling.yaml` | Carbohydrate-active enzyme functions | KEGG, Pfam, CAZy |
| `nitrogen_cycling.yaml` | Nitrogen transformation functions | KEGG, Pfam |
| `antibiotics.yaml` | Antibiotic resistance and biosynthesis | KEGG, Pfam |

## Adding a custom category

1. Copy an existing file: `cp carbon_cycling.yaml my_category.yaml`
2. Edit `name`, `description`, `citations`, and patterns
3. Pass the directory to the tool: `--categories-dir categories/`

The tool loads all `.yaml` files in the specified directory at runtime.
Categories with empty `cazy_subcategories: {}` skip CAZy family matching but still match KEGG and Pfam.
