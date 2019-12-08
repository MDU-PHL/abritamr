# abriTAMR

[![CircleCI](https://circleci.com/gh/MDU-PHL/abritamr.svg?style=svg&circle-token=a54d59b013a30a507621695e738f0a72e47d6969)](https://circleci.com/gh/MDU-PHL/abritamr)

abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on a list of given isolates and collates the results into a table, separating genes identified into functionally relevant groups.


## Command-line tool

`abritamr`:

```bash
abritamr --help
usage: abritamr [-h] [--mduqc] [--qc QC] [--positive_control POSITIVE_CONTROL]
                [--Singularity] [--singularity_path SINGULARITY_PATH]
                [--contigs CONTIGS] [--amrfinder_output AMRFINDER_OUTPUT]
                [--workdir WORKDIR] [--resources RESOURCES]
                [--drug_classes DRUG_CLASSES] [--jobs JOBS] [--keep]

MDU AMR gene detection pipeline

optional arguments:
  -h, --help            show this help message and exit
  --mduqc, -m           Set if running on MDU QC data. If set please provide
                        the MDU QC .csv as further input and an additional
                        output suitable for lims input will be provided.
                        (default: False)
  --qc QC, -q QC        Name of checked MDU QC file. (default:
                        mdu_qc_checked.csv)
  --positive_control POSITIVE_CONTROL, -p POSITIVE_CONTROL
                        Path to positive control sequence - must be set if
                        using -m (default: )
  --Singularity, -S     If using singularity container for AMRfinderplus
                        (default: False)
  --singularity_path SINGULARITY_PATH, -s SINGULARITY_PATH
                        Path to the singularity container for AMRfinderplus,
                        if empty will defualt to shub://phgenomics-
                        singularity/amrfinderplus (default: )
  --contigs CONTIGS, -c CONTIGS
                        Tab-delimited file with sample ID as column 1 and path
                        to assemblies as column 2 (default: )
  --amrfinder_output AMRFINDER_OUTPUT, -o AMRFINDER_OUTPUT
                        Tab-delimited file with sample ID as column 1 and path
                        to amr_finder output files as column 2 (default: )
  --workdir WORKDIR, -w WORKDIR
                        Working directory, default is current directory
                        (default: /home/khhor/MDU/JOBS/ad_hoc/abritamr_in_para
                        llel/20191206)
  --resources RESOURCES, -r RESOURCES
                        Directory where templates are stored (default: /home/k
                        hhor/miniconda3/envs/abritamr/lib/python3.7/site-
                        packages/abritamr)
  --drug_classes DRUG_CLASSES, -d DRUG_CLASSES
                        Path to file (default: /home/khhor/miniconda3/envs/abr
                        itamr/lib/python3.7/site-
                        packages/abritamr/db/refgenes.csv)
  --jobs JOBS, -j JOBS  Number of AMR finder jobs to run in parallel.
                        (default: 16)
  --keep, -k            If you would like to keep intermediate files and
                        snakemake log. Default is to remove (default: False)
```