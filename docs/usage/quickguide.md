# Quick start



`bohra` is a flexible pipeline and allows users to customise the workflows used. This page will provide you a quick start guide to how to trigger the more common features. You can find additional details and options for these and other pipelines [here](../pipelines/overview.md)

If you just want to see some frequently asked questions and commands check out the [cookbook](../cookbook/how_to.md)

Users have the option to supply paired-end fastq (support for ONT coming soon) and/or _de novo_ assemblies as inputs into the pipeline. 

- Paired-end fastq only - where you only have paired-end fastq `bohra` will generate _de novo_ assemblies if required.

- _de novo_ assemblies only - where you only have _de novo_ assemblies `bohra` will not use a reference-based approach for comparative analysis. However, reference-free comparative tools, `ska2` and `mash` are available.

- Both paired-end fastq and _de novo_ assemblies - in situations where you have _de novo_ assemblies already generate, you can supply both sequence types to `bohra`. It will use the supplied _de novo_ assemblies for any steps which require them, potentially saving you time.

Additionally, you can also supply the species value and any optional sample metadata that may be useful.



If you are moving across from bohra version 2 check out the [migration information](../usage/migration.md) 

## 1. Generate your input file

`bohra` requires a single tab-delimited file as input. This can be generated with:

```
bohra generate-input --reads /path/to/reads --contigs /path/to/contigs --outname my_data.txt
```
This will generate a file called `my_data.txt` (defaults to `bohra_input.tsv`) which you can use as input into `bohra`.

**Please note**
You must supply at least one of `--reads` OR `--contigs`. You can supply both - but please be aware that your assemblies must be named in line with your sample/isolate names. For example `sample_name_A.fa`, `sample_name_B.fa` and so on.

Details of input file can be found [here](../usage/inputs.md)

## 2. Choose a pipeline
All pipelines will output a directory with a html summary file and line list results from all sequences in your analysis. This folder is called `report` by default but can be set using the `--report_outdir` flag.

### For quality control or rapid assessment of a dataset

```
bohra run preview -i input_file.tsv -j my_basic_pipeline --cpus N --report_outdir your_report
```

This command will run basic read assessment and also run `mash` to allow you to assess your dataset for poor quality and/or identify outliers which are not suitable for a comparative analysis. 

It can also give you good overview for understanding your pathogen population and determining next steps.

### For comparative analysis (AKA trees/snps)

**`snippy`**

`snippy` is the default comparative analysis tool. 


```
bohra run comparative -i input_file.tsv -j my_snippy_pipeline --cpus N --report_outdir your_report -ref your_reference.fa
```

`bohra` will also accept a gbk reference genome. However, correct `snippy` functioning requires a fasta formatted genome. `bohra` will check the format and generate a `snippy` friendly reference to align to. 

**`ska2`**

This is an alternative to reference based comparative analysis. You can use ska2 in `bohra` if you have assemblies rather than reads, or if you have a mix. Or if you simply want to use a reference free approach.

Please note that when using large datasets, ska2 may take a long time to run. It is advisable to run ska2 on datasets that you already know are likely related to reduce run time and also increase the quality of the comparisons.


```
bohra run comparative --comparative_tool ska -i input_file.tsv -j my_ska_pipeline --cpus N --report_outdir your_report
```

### For AMR and serotyping

`bohra` has serotyping options for many species, the full list can be found [here]. In addition, where available species-specific AMR gene detection and phenotypic predictions may be undertaken.  There is no need for you to specify which typer to use, this will be determined based on species detected in the sequences supplied.


```
bohra run amr_typing -i input_file.tsv -j my_typing_pipeline --cpus N --report_outdir your_report
```

If you have supplied paired-end reads as your inputs, `bohra` will generate assemblies for use in typing and amr processes. The default assembler is `shovill` with `spades`, but users can change this if required. See [advanced usage]().

### Complete pipeline

```
bohra run full -i input_file.tsv -j my_full_pipeline --cpus N --report_outdir your_report
```
This pipeline will run basic sequence assessment, comparative analysis, typing and AMR gene detection.

## 3. Interpretation of report.

The `bohra` pipeline generates a folder for each sequence in the analysis with raw results. In addtion a summary folder with all the combined results (name set with `--report_outdir` or `report` by default) from all sequences in the analysis as well as a html report file, that can be shared.

```
|--seqid_1
      | all results files
|--seqid_2
      | all results file
|--report
      | summarised results files
      | bohra.html
```

- Summary tab (`summary.tsv`) is a collection of key results from the analysis and will include basic sequence metrics, species information and provides information about the quality of the sequence. 
- Species tab (`species.txt`) is a summary of the kraken2 results, indicating the top 3 species identified in the sequences as well as the amount of unclassified. Where you have paired-end and assembly data, the species from both of these sequence types will be provided for comparison.
- Resistome tab (`resistome.txt`) and Virulence (`virulence.txt`) are the raw results of `abritamr`, with genes detected from assemblies grouped by the drug class to which they are assigned.
- Reportable AMR tab (`reportable_amr_mechanisms.txt`) displays the AMR mechanisms that were identified and are classified as relevant for reporting in the context of clinical and/or public health based on the species that is detected in the sequence.
- Core genome stats tab (`core_genome_stats.txt`) details the quality of each alignment. Graphical depiction of the variants across sites in the genome as well as the distribution of alignment metrics are accompanied by a table detailing whether the % alignment. Note that sequences with an alignment % < 2SD from median will be excluded from analysis by default. The behaviour can be changed by using the `--ignore-warnings` flag when running a comparative analysis.

For detailed information about other output files and tables please check [here](../usage/report.md).


## Important considerations.
**Species_expected column**

If you do not require speciation as part of the pipeline and already know the species, you can provide it here. Please note if no speciation is undertaken, `bohra` will use this value to undertake typing and AMR mechanisms/inferrence. If the species in this column is NOT accurate - unexpected results may occur. 

Please note that species detected from sequence trumps the Species_expected value. Were you supply a value in the Species_expected column and also run speciation, the species detected from the sequences will be used to determine the AMR and also serotyping outputs.  If you would like to run bohra and force a particular species, then use the Species_expected column and `--speciation none`.


**Controls**

If you include controls in an analysis you can add an `is_control` column to your input file and supply the type of control. Note that anything designated as a control will have single-sample analyses done (were possible) but will NOT be included in any comparative analyses.


**Tree Annotation and additional metadata RETURNING SOON**

You may also provide additional columns of relevant metadata in your input file. `bohra` will NOT do any data validation on these columns - that is up to the user. But any additional metadata provided will be visible on the tree provided in the report file and the summary table.


**A note on databases**

Many bioinformatics tools require the use of a database or data collection. Where possible and appropriate, `bohra` utilises the databases that come packaged with the tools being used in order to ensure expected behaviour and consistency. However, there are cases were the user will need to supply a database path. 

* `kraken2` - these databases are very large and not easily packaged for distribution with software. So the user will either need to have existing databases available or download them as part of the setup of `bohra`. This information is [here](../installation.md)

* `mlst` comes packaged with a collection of profiles and is ready to use. However, due to changes in licensing, the most up to date profiles cannot be included. As such if you have available to you a current mlst database that is configured for use with `mlst` you can set the `BOHRA_PUBMLST_DB` and `BOHRA_BLAST_DB` environment variables as part of the setup of `bohra`, as described [here](../installation.md)
