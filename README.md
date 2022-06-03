# abriTAMR

[![CircleCI](https://circleci.com/gh/MDU-PHL/abritamr.svg?style=svg&circle-token=a54d59b013a30a507621695e738f0a72e47d6969)](https://circleci.com/gh/MDU-PHL/abritamr)

**_Taming the AMR beast_**

abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on a single (or  list ) of given isolates and collates the results into a table, separating genes identified into functionally relevant groups.

_abriTAMR is accredited by NATA for use in reporting the presence of reportable AMR genes in Victoria Australia._

## New look in v 1.0.0

* Acquired resistance mechanims in the form of point mutations (restricted to subset of species)
* Streamlined output.
* Presence of virulence factors

## Install


abriTAMR requires [AMRFinder Plus](https://github.com/ncbi/amr), this can be installed with `conda`.

abriTAMR comes packaged with an AMRFinder DB consistent with current NATA accreditation. If you would like to use another DB please download it and use the `-d` flag to point to your database.

```
conda create -n abritamr -c bioconda ncbi-amrfinder
conda activate abritamr
pip3 install abritamr
```


## Command-line tool

```bash
abritamr run --help

optional arguments:
  -h, --help            show this help message and exit
  --contigs CONTIGS, -c CONTIGS
                        Tab-delimited file with sample ID as column 1 and path
                        to assemblies as column 2 OR path to a contig file
                        (used if only doing a single sample - should provide
                        value for -pfx). (default: )
  --prefix PREFIX, -px PREFIX
                        If running on a single sample, please provide a prefix
                        for output directory (default: abritamr)
  --species {Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes,Vibrio_cholerae}, -sp {Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes,Vibrio_cholerae}
                        Set if you would like to use point mutations, please
                        provide a valid species. (default: )
  --jobs JOBS, -j JOBS  Number of AMR finder jobs to run in parallel.
                        (default: 16)
```

You can also run abriTAMR in `mdu` mode, this will output a spreadsheet which is based on reportable/not-reportable requirements in Victoria. You will need to supply a quality control file (comma separated) (`-q`), with the following columns:

* ISOLATE
* SPECIES_EXP (the species that was expected)
* SPECIES_OBS (the species that was observed during the quality control analysis)
* TEST_QC (PASS or FAIL)

`--sop` refers to the type of collation and reporting pipeline
* general
  * standard reporting structure for aquired genes, output as reportable and non-reportable
* salmonella
  * Inferred AST based on validation undertaken at MDU

```abritamr mdu --help            
  --qc QC, -q QC        Name of checked MDU QC file. (default: )
  --runid RUNID, -r RUNID
                        MDU RunID (default: Run ID)
  --matches MATCHES, -m MATCHES
                        Path to matches, concatentated output of abritamr
                        (default: summary_matches.txt)
  --partials PARTIALS, -p PARTIALS
                        Path to partial matches, concatentated output of
                        abritamr (default: summary_partials.txt)
  --sop {general,salmonella}
                        The MDU pipeline for reporting results. (default:
                        general)

```
