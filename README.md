# abriTAMR

[![CircleCI](https://circleci.com/gh/MDU-PHL/abritamr.svg?style=svg&circle-token=a54d59b013a30a507621695e738f0a72e47d6969)](https://circleci.com/gh/MDU-PHL/abritamr)q

abriTAMR runs AMRFinderPlus on a list of given isolates and collates the results into a table, separating genes identified into functionally relevant groups.


## Command-line tool

`mdu-amr-detection`:

```bash
mdu-amr-detection --help
usage: mdu-amr-detection [-h] [--mduqc] [--contigs CONTIGS]
                         [--amrfinder_output AMRFINDER_OUTPUT]
                         [--workdir WORKDIR] [--resources RESOURCES]
                         [--drug_classes DRUG_CLASSES] [--jobs JOBS] [--keep]

MDU AMR

optional arguments:
  -h, --help            show this help message and exit
  --mduqc, -m           Set if running on MDU QC data. If set please provide
                        the MDU QC .csv as further input and an additional
                        output suitable for lims input will be provided.
                        (default: False)
  --contigs CONTIGS, -c CONTIGS
                        Tab-delimited file with sample ID as column 1 and path
                        to assemblies as column 2 (default: )
  --amrfinder_output AMRFINDER_OUTPUT, -o AMRFINDER_OUTPUT
                        Tab-delimited file with sample ID as column 1 and path
                        to amr_finder output files as column 2 (default: )
  --workdir WORKDIR, -w WORKDIR
                        Working directory, default is current directory
                        (default: /Users/andersgoncalves/OneDrive - The
                        University of Melbourne/dev/abritamr)
  --resources RESOURCES, -r RESOURCES
                        Directory where templates are stored (default:
                        /Users/andersgoncalves/OneDrive - The University of
                        Melbourne/dev/abritamr/abritamr)
  --drug_classes DRUG_CLASSES, -d DRUG_CLASSES
                        Path to file (default: /Users/andersgoncalves/OneDrive
                        - The University of
                        Melbourne/dev/abritamr/abritamr/db/refgenes.csv)
  --jobs JOBS, -j JOBS  Number of AMR finder jobs to run in parallel.
                        (default: 16)
  --keep, -k            If you would like to keep intermediate files and
                        snakemake log. Default is to remove (default: False)
```