# `linelist`

`linelist` creates a sample-level summary of mechanisms and AMR types from
`amr_status` results.

```sh
abritamr linelist --amr isolate-1.status.csv > isolate-1.linelist.csv
```

`--amr` (`-a`) accepts an `amr_status` result file and defaults to standard
input. The default `--viewtype compact` reports detected mechanisms; choose
`--viewtype full` to include the available drug classes. Use `--genesonly` to
exclude SNPs and `--min-identity` or `--min-coverage` to adjust detection
thresholds. Output is CSV by default; `--format tab` selects tab-delimited
output.
