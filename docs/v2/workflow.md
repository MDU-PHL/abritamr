# Full workflow

## Run each command separately

The individual commands expose each intermediate result. This example runs a
single assembly through detection, status assignment, line listing, matrix
generation, and inference:

```sh
abritamr scan --contigs assembly.fasta --sample-id isolate-1 \
  --species "Salmonella enterica" > isolate-1.scan.csv
abritamr amr_status --amr isolate-1.scan.csv > isolate-1.status.csv
abritamr linelist --amr isolate-1.status.csv > isolate-1.linelist.csv
abritamr matrix --amr isolate-1.status.csv > isolate-1.matrix.csv
abritamr infer --amr isolate-1.scan.csv \
  --reference-folder /path/to/reference-data \
  > isolate-1.infer.csv
```

The status and matrix steps use the default installed reference catalog unless
`--reference-catalog` is provided; `linelist` uses the status results directly.
`infer` requires species rules in its reference folder; inference output may be
empty if no rules are available for the species. To run the entire sequence in
one command, use [`abritamr run`](#run-the-complete-pipeline).

## Run the complete pipeline

For one assembly, provide its identifier and an output directory:

```sh
mkdir -p results
abritamr run --contigs assembly.fasta --sample-id isolate-1 \
  --species "Salmonella enterica" --workdir results --outdir isolate-1
```

This runs scan, status assignment, linelist, matrix, and inference. The output
directory contains `abritamr_scan.csv`, `abritamr_typed.csv`,
`abritamr_linelist.csv`, `abritamr_matrix.csv`, and `abritamr_infer.csv` when
each result is available. `--format tab` writes tab-delimited `.txt` files.

For multiple assemblies, create a tab-delimited file with a `sample_id` column
and a `contigs` column containing assembly paths:

```text
sample_id	contigs
isolate-1	/data/isolate-1.fasta
isolate-2	/data/isolate-2.fasta
```

Then run:

```sh
mkdir -p results
abritamr run --multi samples.tsv --workdir results
```

Each sample gets its own output directory below `results`. Run
`abritamr run --help` to see all workflow options. The `run` command also
accepts existing AMRFinderPlus output with `--amrfinderplus`.
