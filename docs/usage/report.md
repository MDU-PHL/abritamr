# Navigating the `bohra` report

The `bohra` pipeline generates a html report that allows you to share your results quickly and also to interrogate and visualise your analysis in a single unified space.

This report will be different based on the pipeline you ran, with not all result tables and figures visible for all pipelines.

An example report can be found [here](../images/example_report.html)

## Header

The report header provides a summary of the analysis, list user name, pipeline that was run, data, number of sequences in the analysis and the `bohra` version.

## Common tabs

The following analysis results will be available in all pipelines

### Summary

This table presents highlights from your analysis and content will vary depending on the pipeline you ran. It is designed to present the user with information regarding the dataset and highlight outliers in the data that should be removed or further investigated.

**Quality assessment**
The `QA` and `Comment` column represent the results of `bohra` assessment of the sequences. This information should be used as a guide. Where a sequence is green, this indicates that 
- Avg qscore > 30
- Estimated genome depth of coverage (reads only) > 40X
- Number of contigs (assembly only) within 2SD of median number of contigs
- The genome size estimated is within the expected range of the species detected.
- Where alignment has been performed, the % alignment is within 2SD of median % alignment.

Where any of these conditions are not met an orange indicator will draw your attention, with the `Comment` column providing a brief description of what to look for.

Additionally there is a graph button which displays graphical representation of key sequence metrics.

### Sequence assessment

Reads and assembly sequence metrics will be supplied (depending on the input supplied by the user). These tables will report all the key sequence metrics for that input type.

### Speciation

This tab will detail the amount of unclassified reads in each sequence and the top 3 species identified for each sequence using the kraken2 output. Were reads and assemblies are used this data will be lined up together. Addtionally there is a `Species` column which is the `bohra` determined value that is used as input for AMR and serotyping.

### Versions

This table reflects all the tools that were used in the analysis, with references and version numbers.

## AMR and typing tabs

### MLST

This tab outlines the MLST results for all sequences in the dataset.

### Typing

`bohra` undertakes species specific serotyping for 
- _S. enterica_ 
- _L. monocytogenes_ 
- _N. meningiditis_ 
- _N. gonorrhea_ 
- _E. coli_
- _Shigella_ species
- _S. pyogenes_
- _S. aureus_
- _C. auris_
- _Klebsiella_ species

Each serotyper will have its own menu tab (where more than one typer is run). Key columns will be displayed, with raw results available in each sequence directory.

### Plasmid
This is a nested table for each sequence with details of any potential plasmids identified, the species of origin of the plasmid and any AMR genes that were identified.

### Resistome
This table is a combination of the summary output of `abritAMR` and `mob_suite`, identifying AMR genes and SNPs and annotating possible association with a plasmid.

### Reportable_AMR
This table is generated from `abritAMR` and represents the genes and SNPs detected that are clincally relevant/reportable for the species detected.

### Virulence
Similar to the resistome table, this table represents virulence factors identified from `abritAMR`.

## Comparative tabs

### Tree
Were a comparative or preview pipeline was run a tree will be generate and displayed. This tree can be navigate via right clicking on internal nodes to step into the tree. 

### Reference based comparison

Were `snippy` was used as the comparative tool (default) a `Core genome stats` tab will be visible. This table will detail the alignment statisitics and also has a graph for users to identify outliers. Were a sequence is < 2SD from median, it will be annotated in this table. Additionally, a graph with the distribution of variants across the reference genome will also be visible.

### Reference free comparison

Were `ska2` is used as the comparative tool, there will be a `SKA summary stats` file which will provide information about the quality of each pairwise comparison made.

### Distances

The distance matrix is accompanied by two visualisations, a heat map and a graph of the pairwise distance distribution.

### Clusters

`bohra` will also perform heirachical clustering (using single-linkage) based on user supplied thresholds. The table will provide a view of the distance matrices for each cluster and also graphical representation of the inter and intra cluster distances at each threshold.

### Pangenome

Where the `full` pipeline is run, `bohra` generates some pangenome visualisations, indicating the breakdown of core and pan genes and also were clustering is done, a graph showing how the pangenome for each group relates.