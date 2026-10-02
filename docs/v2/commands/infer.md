# `infer`

`infer` applies species-specific genotypic drug susceptibility rules to scan
results. Inference is available only when rules exist for the sample's species.

```sh
abritamr infer --amr isolate-1.scan.csv \
  --reference-folder /path/to/reference-data > isolate-1.infer.csv
```

`--amr` (`-a`) accepts scan results and defaults to standard input. The sample's
species must be present in those results (set it during `scan`).
`--reference-folder` points to the directory containing files named
`02_abritamr_<species>_rules.csv`. Results are CSV by default;
`--reporttype wide` produces a wide rather than the default long format.
`--dflt_result` sets the result for drugs for which no rule matches. Use
`--format tab` for tab-delimited output.

See [Catalog and rules](../catalog-and-rules.md) to create a compatible rules
directory.
