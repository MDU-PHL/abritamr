# abriTAMR

[![CircleCI](https://circleci.com/gh/MDU-PHL/abritamr.svg?style=svg&circle-token=a54d59b013a30a507621695e738f0a72e47d6969)](https://circleci.com/gh/MDU-PHL/abritamr)

**_Taming the AMR beast_**

abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on a single (or  list ) of given isolates and collates the results into a table, separating genes identified into functionally relevant groups.

_abriTAMR is accredited by NATA for use in reporting the presence of reportable AMR genes in Victoria Australia._

* Acquired resistance mechanims in the form of point mutations (restricted to subset of species)
* Streamlined output.
* Presence of virulence factors

## Install

abriTAMR requires [AMRFinder Plus](https://github.com/ncbi/amr), this can be installed with `conda`.

abriTAMR comes packaged with an AMRFinder DB consistent with current NATA accreditation. If you would like to use another DB please download it and use the `-d` flag to point to your database.

```
conda create -n abritamr -c bioconda abritamr
conda activate abritamr
```


## Command-line tool

```bash
abritamr run --help


optional arguments:
  -h, --help            show this help message and exit
  --contigs CONTIGS, -c CONTIGS
                        Tab-delimited file with sample ID as column 1 and path to assemblies as column 2 OR path to a contig
                        file (used if only doing a single sample - should provide value for -pfx). (default: )
  --prefix PREFIX, -px PREFIX
                        If running on a single sample, please provide a prefix for output directory (default: abritamr)
  --jobs JOBS, -j JOBS  Number of AMR finder jobs to run in parallel. (default: 16)
  --identity IDENTITY, -i IDENTITY
                        Set the minimum identity of matches with amrfinder (0 - 1.0). Defaults to amrfinder preset, which is 0.9
                        unless a curated threshold is present for the gene. (default: )
  --amrfinder_db AMRFINDER_DB, -d AMRFINDER_DB
                        Path to amrfinder DB to use (default:
                        /home/khhor/dev/abritamr/abritamr/db/amrfinderplus/data/2021-09-30.1)
  --species {Neisseria,Clostridioides_difficile,Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes}, -sp {Neisseria,Clostridioides_difficile,Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes}
                        Set if you would like to use point mutations, please provide a valid species. (default: )
```

You can also run abriTAMR in `rpeort` mode, this will output a spreadsheet which is based on reportable/not-reportable requirements in Victoria. You will need to supply a quality control file (comma separated) (`-q`), with the following columns:

* ISOLATE
* SPECIES_EXP (the species that was expected)
* SPECIES_OBS (the species that was observed during the quality control analysis)
* TEST_QC (PASS or FAIL)

`--sop` refers to the type of collation and reporting pipeline
* general
  * standard reporting structure for aquired genes, output as reportable and non-reportable
* plus
  * Inferred AST based on validation undertaken at MDU

```
abritamr report --help
 
optional arguments:
  -h, --help            show this help message and exit
  --qc QC, -q QC        Name of checked MDU QC file. (default: )
  --runid RUNID, -r RUNID
                        MDU RunID (default: Run ID)
  --matches MATCHES, -m MATCHES
                        Path to matches, concatentated output of abritamr (default: summary_matches.txt)
  --partials PARTIALS, -p PARTIALS
                        Path to partial matches, concatentated output of abritamr (default: summary_partials.txt)
  --sop {general,plus}  The MDU pipeline for reporting results. (default: general)
```
