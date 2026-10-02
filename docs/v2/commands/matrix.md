# `matrix`

`matrix` creates a sample-level presence/absence table from `amr_status` results.

```sh
abritamr matrix --amr isolate-1.status.csv > isolate-1.matrix.csv
```

`--amr` (`-a`) accepts an `amr_status` result file and defaults to standard
input. The default `--facet abritamr_subclass` groups results by AMR subclass;
choose `abritamr_class` or `gene` to group by class or gene. The selected
`--reference-catalog` determines the class and subclass columns. Use
`--min-identity` and `--min-coverage` to adjust detection thresholds, and
`--format tab` for tab-delimited output.
