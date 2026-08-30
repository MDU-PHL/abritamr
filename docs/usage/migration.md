# Migration from bohra v2 to bohra v3


## Summary of key changes
For users who are were using previous versions of bohra major usage changes are 
1. [New input file structure](../usage/overview.md#input-file), which allows for addition of user supplied metadata, species expected as well as reads and assemblies. It is possible to convert existing version 2 input files ([see below](#convert-bohra-v2-input-file-to-bohra-v3))
2. Each pipeline has its own command now

```
bohra run --help

amr_typing   Help for the amr_typing pipeline.
assemble     Help for the assemble pipeline.
basic        Help for the basic pipeline.
comparative  Help for the comparative pipeline.
full         Help for the full pipeline.
tb           Help for the tb pipeline.
```

2. When undertaking a comparative analysis bohra v3 will do clustering using heirarchical clustering (average, complete or single linkage) with user defined SNP thresholds.

3. Support for using assemblies alone as input (without having to have fastq files) is also a new feature of bohra v3. You can now supply assemblies (from short-read Illumina or even ONT) and use [ska](https://github.com/bacpop/ska.rust/) to maeasure distances and build trees. This can make comparative analysis quicker and easier.

4. You can also supply a mixture of reads and assemblies when using ska or mash (_possible but interpret these results with caution_).

5. Sample associated metadata is now supported. Where you have metadata, such as source, geography etc; you can supply this information in the input file and it will be presented in result tables and used to annotate the tree.

6. Additional _in silico_ serotyping is also available, ShigaPass and sonneitype are now run where _Shigella_ species is detected.

7. New look report html with some additional visualisations and more information about what was run. EXAMPLES COMING SOON

## Convert bohra v2 input file to bohra v3

bohra v3 uses a single input file, rather than multiple that were previously used. If you have existing bohra version to input files you can convert them with csvtk (bohra command will shortly be available).

1. You have an input file with only reads

```
conda activate bohra
csvtk -t add-header -n 'Isolate,r1,r2' reads.tab > new_bohra_input.txt
```

2. If you have both reads and contigs inputs for bohra version 2. 

It is a good idea if you have both reads and contigs to supply the contigs in the input file - this will prevent the time consuming step of assemblng genomes.

```
conda activate bohra
csvtk -t -H join --outer-join -f1 reads.tab contigs.tab | csvtk -t add-header -n 'Isolate,r1,r2,assembly' > new_bohra_input.txt
```
Please note when running the above command **the order you supply the files is important**. It is recommended that you check the resulting file to make sure that this has worked as expected

```
csvtk -t pretty new_bohra_input.txt | less -S
```
Should result in 

| Isolate | r1 | r2 | assembly |
|:--- | :--- | :---| :--- |
| seq1 | /path/to/seq1_read1.fastq.gz| /path/to/seq1_read2.fastq.gz | /path/to/seq1_contig.fa|
| seqn | /path/to/seqn_read1.fastq.gz| /path/to/seqn_read2.fastq.gz | /path/to/seqn_contig.fa|


3. Create an input file from just contigs

```
conda activate bohra
csvtk -t add-header -n 'Isolate,assembly' contigs.tab > new_bohra_input.txt
```
If you would like to generate a completely new input file check out the instructions [here](../usage/overview.md#input-file).

## How to run a pipeline

bohra version 2 pipelines can be run in bohra v3 

| bohra v2 command | bohra v3 command | bohra v3 input |
|:--- | :---| :--- | 
| `bohra run -p full` | `bohra run full` | paired-end fastq and/or assemblies |
| `bohra run -p snps` | `bohra run comparative`| paired-end fastq and/or assemblies |
|`bohra run -p preview` | `bohra run comparative --comparative_tool mash` | paired-end fastq and/or assemblies |
|`bohra run -p amr_typing` | `bohra run amr_typing`| paired-end fastq and/or assemblies |
| Not available | `bohra run tb` | paired-end fastq | 



## Running with just assemlbies

If you would like to run bohra v3 with just assemblies

1. Create an input file with only assembly paths.
2. Choose a pipeline to run:
    - `amr_typing`
    - `full`
    - `comparative` (this will only run comparative tools - there will be no typing or AMR)
**amr_typing**
```
bohra run amr_typing -i new_bohra_input.txt --cpus X
```
(don't forget to add `--kraken2_db` if you have not got `KRAKEN2_DEFAULT_DB` set)

**full (or comparative)**
```
bohra run full (or comparative) -i new_bohra_input.txt --comparative_tool ska (or mash)
```

