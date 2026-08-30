# Troubleshooting and FAQ

If you run into an error or other issue that you are unsure of or is not clear - please create an [issue](https://github.com/MDU-PHL/bohra/issues). If we are not aware of a problem it won't get fixed.


You can also check out the advanced options [here](../usage/advanced.md).


## I only have reads - can I still use `bohra`?

Yes you can. Primarily most of the processes in `bohra` are based on pe-reads and were assemblies are required they will be generated and stored in the sample directory for future use.

## What if I only have assemblies - will `bohra` still run?

Yes it will. However, it is important to be aware that any process which require reads, such as some serotyping, snippy and some types of AMR mechanisms may not be run if there are no reads. If you would like to do a comparative analysis and only have assemblies, you can use `--comparative_tool mash` or `--comparative_tool ska`.

## Why does `bohra` run kraken on reads and assemblies?

For quality assessment of your data, if the most abundant species in the pe-reads is different to the most abundant species detected in the assembled genomes, this indicates that there is something strange in your sequence data. It could be contamination, it could be low sequence depth and it is important to further investigate these cases.



## I have an error `Output directory already exists`

`bohra` does not blindly overwrite existing aggregated and summarised report folder. 
- If you DO want to override previous report folder you can follow the suggested fix of adding `--replace-report` to your command.
- If you DO NOT want to override the previous report folder add `-ro YOUR_NEW_REPORT_NAME` to your command

## `bohra` is not including some of my sequences in the core genome analysis

To prevent pipeline failure and inclusion of strange sequences in a core genome analaysis, `bohra` will remove any sequences that have a % alignment to the reference < 2SD from median. It will provide a comment stating that these sequences are outside the expected range. 
- If you want to include all the sequences you can add `--ignore_warnings` to your command. The warning/information will still appear in the tables, but the sequence will NOT be filtered out of your core genome analysis.

## What are the quality thresholds used in `bohra`?

The quality assessment used in `bohra` are largely provided as a guide. The only metrics which will directly affect your analysis are 
- if the filesize < 2MB only sequence assessment will be done - no other tests/processes will be run (ie no assembly or alignment)
- if the % alignment is < 2SD from median (see [above](../cookbook/faq.md#bohra-is-not-including-some-of-my-sequences-in-the-core-genome-analysis))

Other metrics are simply used to direct the user to investigate the sequence further. `bohra` does not make any assumptions about the suitability for inclusion in an analysis. This is up to the user. The only criteria that `bohra` imposes are to prevent pipeline failure.

| Metric | Threshold | Input type | 
| :---:| :---:| :---:| 
| Filesize | > 2000000 MB | pe-reads|
| Genome size | Genome size as assessed from PE-reads and/or number of base pairs in assembly within range for expected species | pe-reads and assemblies
| Estimated depth of coverage | > 40x | pe-reads |
| % alignment | < 2SD median dataset alignment | pe-reads with `snippy`|
| number contigs | +/- 2SD median | assembly|

## How is genome size calculated?

Genome size is estimated from both pe-reads as well as assemblies.

- genome size from assemblies is simply the number of bp in the assembly
- genome size from pe-reads is calculated using `kmc`

## Where does GC content calculated from?

Primarily this is provided from pe-reads. However, if you only supply assemblies, this will be used to present the GC content.