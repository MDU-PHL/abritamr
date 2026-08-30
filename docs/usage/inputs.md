# Inputs

## Input file

`bohra` requires a tab-delimited file detailed below.


| Column name | Description | Required? | 
| :---: | :---: | :---: | 
| Isolate | This is the name of the sequence or sample and will appear throughout `bohra` outputs. It must be unique.| Yes|
| r1 | The path to read 1 | If an assembly file is not supplied you must supply reads | 
| r2 | The path to read 2 | If an assembly file is not supplied you must supply reads |
| assembly | The path to the assembly for the isolate | If reads are not supplied you must supply an assembly file |
| Species_expected | The expected species of the sample or 'control'. the proper species name (not the _ joined name from amrfinder) | No |


**`bohra` can generate the input file for you**

If you have 


* Paths to your reads and/or contigs

* (Optional) A table with a list of isolates and other data (species or other metadata) (column 'Isolate' must be included) 

`bohra` can generate the input file for you. 

## Reference file

If you are running `bohra` using one of the comparative pipelines with `snippy` as the comparative tool, then you will need to supply a reference file. `bohra` can accept any format, `.fa*`, `.gbk` and `gz` files. File format will be converted to `.fasta` for use with `snippy`.

