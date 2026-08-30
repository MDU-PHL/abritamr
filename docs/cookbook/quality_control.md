# Sequence quality control

Assessing quality of your sequence data is vitally important to facillitate appropriate interpretations of genomic analysis results.

`bohra` provides some basic quality checks which can be useful to determine the quality of your data (summary tab of the analysis report). The pipeline will NOT remove anything from single sample analysis, but will flag sequences which display features that could detrimentally impact the quality of the analysis. Where a sequence is identified as an outlier in a comparative analysis (core genome analysis) it will be removed to prevent the pipeline from failing. However, you can provied the `--ignore-warnings` flag to prevent this behaviour.

**Summary table**

The `bohra` pipeline outputs a html report with all results for the pipeline that was run. The summary table will be the first or second tab you see when opening up the file in your browser.

## Key metrics

Depending upon your input data and the pipeline run some or all of the following may be available to you.

### GC content

GC content (the proportion of bases that are G or C) can be calculated from fastq or fasta files. It can be a marker of species, for example _Listeria monocytogenes_ often has a low GC content (~37% - 39%), whilst _M. tuberculosis_ can have a higher GC content (~64%). If the GC content of your sequence does not correlate with the species that is expected or detected, it may indicate contamination or the sequence is not what you thought it should be.

### Number of reads

The number of reads presented by the `bohra` pipeline is the estimated number of paired reads. This value are only relevant if they are extremely low or extremely high and may indicate that the quality of the sequencing is poor. The total number of bases in your input is presented as the Yeild of sequencing and by itself is not terribly useful. This metric is used to estimate the average depth of coverage of your sequencing.

### Insert size

The insert size value calculated by `bohra` (calculated for all pipelines that have paired-end fastq and were assemblies are provided or generated) is an estimation of the length of sequence covered by a pair of reads. Ie R1 + sequence + R2. Ideally this value should be > 2x read length, where this value is lower than that, you may see issues in the quality of results generated from assemblies (serotyping and detection of AMR mechanisms).

### IQR depth coverage

The IQR of the depth of coverage calculated by `bohra` indicates the variability in the depth of sequencing coverage. Where the IQR > average depth of coverage, this indicates that there is significant areas of the sequence with 0 depth of coverage. This may cause issues with both assembly based downstream tools as well as alignment based (core genome or Mtb AMR), due to missing data.

### Depth

For paired-end sequences, the depth of coverage is estimated from the estimated genome size and the number of bases in the sequences. This value should be sufficiently high to ensure that there is enough sequencing depth to generate quality _de novo_ assemblies 

### Assembly quality

Assembly quality can be assessed by looking at the N50, number of contigs and the assembly length. The **assembly length** is indicative of the assembled genome size and should reflect the expected genome size of the species present in the sequence. The **number of contigs** is often used as a key metric to assess the quality of a _de novo_ assembly. However, the use of hard and fast thresholds for this metric can be problematic. Generally, a high number of contigs (>1000) can indicate that the assembly quality is poor, but this number can be influenced by the library preparation methods used to generate the fastq, the assembler that is used (some assemblers favour accuracy over contiguity) and species. Furthermore, the **N50** can also be useful to assess that overall assembly quality. This can be used to indicate the contiguity of the sequence, were the N50 is low, this would indicate lots of small fragments, whilst a high number here would be indicative of more contigous sequence. Thus when assessing the quality of an assembly, it is important to consider all  3 of these metrics in context of the dataset you are investigating. 

