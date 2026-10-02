# Species detection

When scanning an assembly, abriTAMR can try to identify its species if you do
not provide one with `--species` (`-sp`). It uses Sourmash to compare the
assembly against the species signatures bundled with abriTAMR. The closest
matching signature is used when the similarity is at least 0.01; otherwise the
species is reported as `unknown`. The current implementation builds the query
from the first sequence record in the assembly, so the result should be treated
as an estimate rather than a definitive identification.

If you know the species, provide it explicitly:

```sh
abritamr scan --contigs assembly.fasta --species "Salmonella enterica"
```

The supplied species is used in preference to automatic identification. This
is especially important when scanning existing AMRFinderPlus output, where
abriTAMR cannot infer a species from an assembly. Species information is used
to select the AMRFinderPlus organism and to support species-specific AMR
annotation and inference. Check the reported species and supply `--species`
when the estimate is missing or incorrect.
