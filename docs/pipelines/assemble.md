# assemble pipeline


This workflow will simply generate assemblies from paired-end fastq, run basic genome annotation with `prokka` and assess the quality of both the input reads and the resulting assemblies. This workflow forms the basis for amr, typing and pangenome analysis.

```mermaid
flowchart LR
fastq --> assembly --> annotation --> sequence_assessment
assembly --> speciation
fastq --> sequence_assessment --> report
fastq --> speciation --> report

```

This pipeline should be used if you would like to generate _de novo_ assemblies from paired-end reads. You can choose from `spades`, `skesa`, `shovill + skesa` or `shovill + spades` (default)

```
bohra run assemble -i input_file.tsv -j my_assembly_pipeline -a shovill_skesa --cpus X
```
