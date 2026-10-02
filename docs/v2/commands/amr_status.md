# `amr_status`

`amr_status` reads `scan` results and adds abriTAMR reporting status, AMR type,
and criteria metadata to each mechanism.

```sh
abritamr amr_status --amr isolate-1.scan.csv \
  --reference-catalog /path/to/01_abritamr_reference_gene_catalog.csv \
  > isolate-1.status.csv
```

The scan result can also be piped on standard input:

```sh
abritamr scan --contigs assembly.fasta --sample-id isolate-1 \
  --species "Salmonella enterica" |
  abritamr amr_status > isolate-1.status.csv
```

`--amr` (`-a`) accepts a scan result file and defaults to standard input.
`--reference-catalog` selects the catalog used for classification. Use
`--format tab` for tab-delimited output. The status result is the input for
[`linelist`](linelist.md) and [`matrix`](matrix.md).
