[![CI](https://github.com/MDU-PHL/abritamr/workflows/CI/badge.svg)](https://github.com/MDU-PHL/abritamr/actions)
[![GitHub release](https://img.shields.io/github/release/MDU-PHL/abritamr.svg)](https://github.com/MDU-PHL/abritamr/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Conda](https://img.shields.io/conda/dn/bioconda/abritamr.svg)](https://anaconda.org/bioconda/abritamr)
[![DOI](https://zenodo.org/badge/209921768.svg)](https://zenodo.org/badge/latestdoi/209921768)

<img src="documentation/abriTAMR_logo.jpg" height="100">

## Taming the AMR beast

_abriTAMR_ is an AMR gene/variant detection 
pipeline  that runs AMRFinderPlus 
on bacteiral genome assemblues,
categorises the variants into 
reportable/not-reportable,
and sorts into clinical drug classes.

_abriTAMR_ is accredited by NATA for use 
in identifying the presence of reportable 
AMR genes the MDU PHL in Victoria, Australia.

## Installation

```
% conda create -n abritamr -c bioconda abritamr
% conda activate abritamr
% abritamr --version

```

## Quick start
```
% abritamr run -c genome.fasta

% ls abritamr
abritamr.log         summary_partials.txt
amrfinder.out        summary_virulence.txt
summary_matches.txt  update_abritamr_db.log

% cat abritamr/abritamr.txt
Isolate	Methicillin	Tetracycline	Tigecycline	Beta-lactam	Penicillin resistance (Staphylococcus aureus)
abritamr	mecA,mecR1^	tet(38)*	mepA*	blaI*,blaR1*	blaZ
```

## Input

The `-c` option is used for input, and 
can accept 2 types of files:

1. a single FASTA file, usually a genome assembly:
```
>contug001
AGTCTCGATATGCTATAGGCTTATATATAT
ATGCTATAGGCTTATATATATTTATATCTT
>contig002
CGATATGCTATAGGCTTATATATATTTATA
...
```

2. a TSV file for multiple FASTA files:
```
ID1 <tab> /path/to/assembly1.fasta
ID2 <tab> /path/to/second/file.fna
...
```

## Running

```
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
                        /<path_to_installation>/abritamr/abritamr/db/amrfinderplus/data/2021-09-30.1)
  --species {Neisseria,Clostridioides_difficile,Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes}, -sp {Neisseria,Clostridioides_difficile,Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes}
                        Set if you would like to use point mutations, please provide a valid species. (default: )
```

You can also run abriTAMR in `report` mode, this will output a spreadsheet which is based on reportable/not-reportable requirements in Victoria. You will need to supply a quality control file (comma separated) (`-q`), with the following columns:

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

## Output

### `abritAMR run` 

Outputs 4 summary files and retains the raw AMRFinderPlus output for each sequence input.

1. `amrfinder.out` raw output from AMRFinder plus (per sequence). For more information please see AMRFinderPlus help [here](https://github.com/ncbi/amr/wiki/Interpreting-results) 

2.  `summary_matches.txt` 
  * Tab-delimited file, with a row per sequence, and columns representing functional drug classes 
  * Only genes recovered from sequence which have >90% coverage of the gene reported and greater than the desired identity threshold (default 90%). 
    
    I. Genes annotated with `*` indicate >90% coverage and > identity threshold < 100% identity.
    
    II. No further annotation indicates that the gene recovered exhibits 100% coverage and 100% identity to a gene in the gene catalog.
    
    III. Point mutations detected (if `--species` supplied) will also be present in this file in the form of `gene_AAchange`.

3. `summary_partials.txt`
  * Tab-delimited file, with a row per sequence, and columns representing functional drug classes 
  * Genes recovered from sequence which have >50% but <90% coverage of the gene reported and greater than the desired identity threshold (default 90%). 

4. `summary_virulence.txt`
  * Tab-delimited file, with a row per sequence, and columns representing AMRFinderPlus virulence gene classification
  * Genes recovered from sequence which have >50% coverage of the gene reported and greater than the desired identity threshold (default 90%). 

      * Genes recovered with >50% but <90% coverage of a gene in the gene catalog will be annotated with `^`.
      * Genes annotated with `*` indicate >90% coverage and > identity threshold < 100% identity.

4. `abritamr.txt`
  * Tab-delimited file, combining `summary_matches.txt`, `summary_partials.txt`, `summary_virulence.txt` with a row per sequence, and columns representing AMRFinderPlus virulence gene classification and/or functional drug classes.
  * Genes recovered from sequence which have >50% coverage of the gene reported and greater than the desired identity threshold (default 90%). 

      * Genes recovered with >50% but <90% coverage of a gene in the gene catalog will be annotated with `^`.
      * Genes annotated with `*` indicate >90% coverage and > identity threshold < 100% identity.

### `abritamr report` 

will output spreadsheets `general_runid.xlsx` (NATA accredited) or `plus_runid.xlsx` (validated - not yet accredited) depending upon the sop chosen.

* `general_rundid.xlsx` has two tabs, one for matches and one for partials (corresponding to genes reported in the `summary_matches.txt` and `summary_partials.txt`). Each tab has 7 columns 

| Column | Interpretation |
|:---: | :---: |
| MDU sample ID | Sample ID |
|Item code | suffix (MDU specific) |
| Resistance genes (alleles) detected | genes detected that are reportable (based on species and drug classification)|
| Resistance genes (alleles) det (non-rpt) | other genes detected that are not not reportable for the species detected.
| Species_obs | Species observed (supplied in input file) |
| Species_exp | Species expected (supplied in input file) |
| db_version | Version of the AMRFinderPlus DB used |

* `plus_runid.xlsx` output is a spreadsheet with the different drug resistance mechanims and the corresponding interpretation (based on validation of genotype and phenotype) for drug-classes relevant to reporting of anti-microbial resistance in _Salmonella enterica_ (other species will be added as validation of genotype vs phenotype is performed).

* Ampicillin
* Cefotaxime (ESBL) 
* Cefotaxime (AmpC)
* Tetracycline
* Gentamicin
* Kanamycin
* Streptomycin
* Sulfathiazole
* Trimethoprim
* Trim-Sulpha
* Chloramphenicol 
* Ciprofloxacin
* Meropenem 
* Azithromycin
* Aminoglycosides (RMT)
* Colistin 

## References

* [AMR FInder Plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)
* [Abricate](https://github.com/tseemann/abricate)

## Feedback

File questions, bugs, or ideas on the 
[Issues page](https://github.com/MDU-PHL/abritamr/issues)

## License

[GPLv3](https://raw.githubusercontent.com/MDU-PHL/abritamr/master/LICENSE)

## Citation

Sherry, N.L., Horan, K.A., ... , Seemann, T. 
_An ISO-certified genomics workflow for identification and surveillance of antimicrobial resistance_
*Nat Commun* 14;60 (2023). 
[DOI:10.1038/s41467-022-35713-4](https://doi.org/10.1038/s41467-022-35713-4)
[PMID:36599823](https://pubmed.ncbi.nlm.nih.gov/36599823/)


## Authors

* Kristy Horan
* [Torsten Seemann](https://tseemann.github.io)
* Norelle Sherry
