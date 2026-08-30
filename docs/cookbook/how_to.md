# How do I

## Generate an input file?
1. If you have only pe-reads
```
bohra generate-input --reads /path/to/reads  --outname my_data.txt
```
2. If you have pe-reads AND pre-made assemblies

```
bohra generate-input --reads /path/to/reads --contigs /path/to/contigs  --outname my_data.txt
```
3. If you have only pre-made assemblies

Please be aware if this is the case only processes which require assemblies will be undertaaken. No pe-read based process will be run
```
bohra generate-input --contigs /path/to/contigs  --outname my_data.txt
```

(you can rename your input as you please - the default input file if no name is supplied is `bohra_input.tsv`)

- You will need a folder with all of your sequences in it.
- `bohra` will match read 1 and read 2 and use the string before the first `_`  as the isolate identifier.
- If your input sequences are assemblies, `bohra` will use the string before the first `.` as the isolate identifier. 
- You do not need to have 
- **Note** if these do not match (ie the name of the assembly is not the same as the reads) `bohra` will potentially create multiple entries.

## Run `mash` to assess dataset?

```
bohra run preview -i input_file.tsv -j my_basic_pipeline --cpus N --report_outdir your_report
```

## Run the full pipeline?

```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa
```

## Change the cluster threshold?
```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa -ct 6,11,21
```
## Change the cluster method
```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa -ct 6,11,21 --cluster_method complete
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa -ct 6,11,21 --cluster_method average
```
There are also `centroid`, `median`, `ward` and `weighted` as possible clustering methods. `single` is the default.

## Use ska instead of `snippy`?
```
bohra run full --comparative_tool ska -i input_file.tsv --cpus N 
```

## Just run amr and serotype?
```
bohra run amr_typing -i input_file.tsv --cpus N
```

## Remove recombination?

```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa --gubbins
```

## Remove (or add a new sequence)?

If you want to remove or add a new sequence to an analysis, simply remove (or add) the new information to the input file and re-rerun your command in the same directory. The pipeline will only re-run process for which there is a new input, there is no need to create a new folder each time. 
- New trees will be run because the input to the tree process will be different
- Sequence assessment and single sample processes (such as assembly, MLST, AMR) will not be re-run on existing samples.

**NOTE** 
In your working directory the pipeline will create a folder called `work`. This is essentially the cache for your all of your analyses. If you remove this folder for any reason - your whole analysis will re-run the next time you trigger it. There will be no history or cache of results. 

## Rerun with different settings?

If you would like to modify a setting, for example modify `--fuzzy_core_prop` setting, you can do this in the same directory as you previously ran in. It is recommended that you do not overwrite the existing report directory, so as to keep record of what your results looked like with different conditions. Similar to adding or removing a sequence, only processes directly impacted by your change will be re-run.