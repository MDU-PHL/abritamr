# The home of abritAMR

**Taming the AMR beast**



_abriTAMR_ is an AMR gene/variant detection pipeline  that runs AMRFinderPlus on bacterial genome assemblies, prioritises the mechanisms based on the species from which the sequences were derived. Additionally, abritAMR will generate genomic DST, based on [AMRrules](https://github.com/AMRverse/AMRrules) or abritAMR (_S. enterica_) for a selection of species.


_abriTAMR_ is accredited by NATA for use in identifying the presence of reportable AMR genes the MDU PHL in Victoria, Australia.

Usage instructions can be found [here](usage/quickguide.md)

## Installation


### Recommended (conda or mamba)


```
% conda create -n abritamr -c bioconda abritamr
% conda activate abritamr
% abritamr --version
```


## Citation

Sherry, N.L., Horan, K.A., ... , Seemann, T. 
_An ISO-certified genomics workflow for identification and surveillance of antimicrobial resistance_
**Nat Commun** 14;60 (2023). 
[DOI:10.1038/s41467-022-35713-4](https://doi.org/10.1038/s41467-022-35713-4)
[PMID:36599823](https://pubmed.ncbi.nlm.nih.gov/36599823/)


## Authors

* Kristy Horan
* [Torsten Seemann](https://tseemann.github.io)
* Norelle Sherry
* CHarlie Higgs (logo design)
