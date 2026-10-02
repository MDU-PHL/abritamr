# Reference catalog and rules

abriTAMR uses a reference gene catalog to annotate AMRFinderPlus hits and
species-specific rule files for genotypic drug susceptibility inference.

## Generate the catalog and rules together

`generate_database` builds both from their configured sources and saves the
files to `abritamr_db` in the current directory by default:

```sh
abritamr utils generate_database --output-dir reference-data
```

The catalog is written as
`reference-data/01_abritamr_reference_gene_catalog.csv`; rules are written as
species-specific files such as `02_abritamr_Salmonella_enterica_rules.csv`.
These files can be used with `--reference-catalog` and `--reference-folder`,
respectively:

```sh
abritamr scan --contigs assembly.fasta \
  --reference-catalog reference-data/01_abritamr_reference_gene_catalog.csv
abritamr infer --amr scan.csv --reference-folder reference-data
```

Generated reference data should be reviewed and kept versioned alongside
analysis results to ensure the same definitions are used when results are
reproduced.

## Generate either component separately

To build only the reference gene catalog:

```sh
abritamr utils generate_gene_catalog --output-dir reference-data
```

To build rules, specify the catalog to use:

```sh
abritamr utils generate_rules \
  --catalog reference-data/01_abritamr_reference_gene_catalog.csv \
  --output-dir reference-data
```

`generate_rules` includes AMRverse rules by default and adds the local
inference definitions. Use `--species` to restrict rule generation to a
configured species, `--evidence_grade` to select the minimum evidence grade,
and `--no-amrrules` to omit AMRverse rules. `--inference-definitions` selects
the local inference-rule CSV; `--output-dir` selects the destination for the
generated rules.

The catalog-generation command accepts `--class-definitions` and
`--amrtyping-definitions` to supply the class curation and reporting criteria
CSV files. All generator options and defaults are available with
`abritamr utils <command> --help`.
