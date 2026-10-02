# `scan`

`scan` runs AMRFinderPlus on assembly input, then annotates the detected
mechanisms with classes from the reference gene catalog. Existing AMRFinderPlus
output can be supplied instead of an assembly.

```sh
abritamr scan --contigs assembly.fasta --sample-id isolate-1 \
  --species "Salmonella enterica" > isolate-1.scan.csv
```

Or provide existing AMRFinderPlus output:

```sh
abritamr scan --amrfinderplus amrfinder.out --sample-id isolate-1 \
  --species "Salmonella enterica" > isolate-1.scan.csv
```

`--contigs` (`-c`) and `--amrfinderplus` (`-afp`) accept one or more files.
`--sample-id` (`-s`) identifies the sample in the results. If omitted, the
input path is used as the sample identifier. `--species` (`-sp`) is optional;
when omitted for an assembly, species is guessed. Species identification is
needed for species-specific SNP detection and inference; supply `--species` if
automatic identification is unavailable or incorrect.

The CSV output contains one row per detected mechanism and includes the sample
ID, species, catalog annotations, and AMRFinderPlus match details. The default
minimum identity is 0.9 and minimum coverage is 0.5; adjust these with
`--min-identity` and `--min-coverage`. Use `--threads` to set the number of
CPUs. `--reference-catalog` selects the annotation catalog and `--species_rules`
selects species-specific AMRFinderPlus rules.

To continue with status assignment, use the resulting file as the input to
[`amr_status`](amr_status.md).
