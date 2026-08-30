# Advanced usage

## Species

`bohra` requires species information to undertake species-specific AMR and serotyping. However, in some situations it may not be feasible or desirable to run speciation within the pipeline. You can provide a pre-determined `Species_expected` value in your input file and using `--speciation none`. Please note though that if the species you supply is NOT what is actually in the sequences you may get unexpected behaviour which can include pipeline failure.

## Ignore the % alignment

If you use the `--ignore_warnings` flag all sequences supplied in your analysis will be included in any comaprative analysis that is done. Please be aware that if your dataset has sequences in it that are not suitable (ie wrong species, or otherwise exhibit poor alignment) you may encounter unexepected behaviour, including pipeline failure. Furthemore any interpretations made were there are outliers may be inappropriate.

## Assembly options

The default assembler for `bohra` is `shovill + spades`. You can select from `shovill + skesa` or `skesa` or `spades`. 

You can also adjust various assembly settings including

- `--min_contig_length` - default is 200 bp, but you can adjust this to your needs
- `--spades_args` - when running `--assembler spades` you can supply strings setting your desired options. For example ` --spades_args '--cov-cutoff auto'`

## Trimming

By default `bohra` does not trim reads. However, if your data is from public sources or older and you are unsure you can add `--trim` to your command and the reads will be trimmed and the trimmed data used as input to all downstream processes. Please note that this will lead to duplication of your reads - so ensure that you have enough storage space available to meet this need.

## MLST 

MLST settings can also be modified, in addition to supplying a `MLST_DBDIR` you can also 
- Exclude schemes using `--mlst_exclude scheme_name_to_exclude`. It is possible to multiples of these `--mlst_exclude scheme_name_to_exclude --mlst_exclude another_scheme_name_to_exclude`

- Output novel allele sequences into the sample directory for the sample they were identified in using `--novel_mlst novel_name` 


## Comparative analysis

### snippy

You can modify  `snippy` arguments as well as define the proportion of core genome to using `core_snp_filter`.

- minimum read mapping quality with `--minmap` (default 60)
- minimum base quality with `--basequal` (default 13)
- minimum QUALITY in VCF column 6 `--minqual` (default 100)
- minimum proportion for variant evidence with `--minfrac` (default AUTO)
- define proportion of core genome to use with `--fuzzy_core_prop` (default 1.0)

You can modify additional `snippy` settings by providing the option and value to `--snippy_args`.

## ska2

You can modify 
- minimum frequency for variant calling with `--ska_minfreq` (default 0.9)
- provide additional arguments for alignment `--ska_alnargs`
- kmer size for ska2 with `--ska2_kszise` (default 31)