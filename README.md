# bwt_filtering_pipeline
The docker image for NGS reads reference alignment pipeline based on [bowtie-tools](https://github.com/ivasilyev/bowtie-tools) repository.

## Input formats
By default, the pipeline consumes two special formats of input files containing tab-delimited text linkages:
- SAMPLEDATA
```
sample_name_11	/path/to/unpaired_sample_file_11
sample_name_12	/path/to/unpaired_sample_file_12
...
sample_name_21	/path/to/paired_sample_file_21.1	/path/to/paired_sample_file_21.2
sample_name_22	/path/to/paired_sample_file_22.1	/path/to/paired_sample_file_22.2
```
- REFDATA
```
/path/to/sequence.fasta	/path/to/colorspace_bowtie_mask	/path/to/regular_bowtie2_mask	/path/to/samtools_index.fai	/path/to/samtools_lengths.genome	/path/to/samtools_annotation.txt
```

The REFDATA may also be generated using the `cook_the_reference.py` script. All required software is included. 

## Main workflow
The pipeline employs `multiprocessing` features through the queueing of multiple samples per time. However, the queueing algorithm is synchronized, so the next step queue shall not be processing while some samples are still running in previous queue. 

Process description
1. Align all samples on the first REFDATA;
2. Align all non-mapped reads left from the first step on the second REFDATA;
3. Convert sequence alignment map files and extract coverage

Due to multithreading implemented in [bowtie](http://bowtie-bio.sourceforge.net)/[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) onboard, the relevant queues are sequental. Other participating applications do not use multiple CPU cores, so the 3rd step queue is parallel.

# metaphlan2_pipeline
WIP

# metaphlan2_pipeline
WIP
