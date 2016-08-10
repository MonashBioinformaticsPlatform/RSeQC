# Release history

## RSeQC v2.6.4

* Two dependency packages bx-python and pysam are **not** shipped with ## RSeQC starting from v2.6.4. 
* Users could install ## RSeQC using pip: **pip install ## RSeQC**. bx-python and pysam will be installed automatically if they haven’t been installed before.

## RSeQC v2.6.3

* Fix a bug in "read_quality.py" that does not return results if input file containing less than 1000 reads.
* ## Remove "## RPKM_count.py" as it generates erroneous results especially for longer reads. Use "FPKM_count.py" instead.
* "bam_stat.py" prints summary statistics to STDOUT. 

## RSeQC v2.6.2

* fix bugs in "insertion_profile.py", "clipping_profile.py", and "inner_distance.py "

## RSeQC v2.6.1

* Fix bug in "junction_annotation.py" in that it would report some "novel splice junctions" that don't exist in the BAM files. This happened when reads were clipped and spliced mapped simultaneously. 
* Add FPKM.py. FPKM.py will report "raw fragment count", "FPM" and "FPKM" for each gene. It does not report exon and intron level count. 

## RSeQC v2.6

Add 3 new modules:

* deletion_profile.py
* insertion_profile.py
* mismatch_profile.py

## RSeQC v2.5

* read_duplication.py:

 * add '-q' option filter alignments with low mapping quality.
 * Fix bug related to the labels of right Y-aixs
 
* bam2fq: add '-c' option to call 'gzip' command to compress output fastq file(s).
* divide_bam.py: add '-s' option, skipped unmapped reads.
* clipping_profile.py:

 * add '-q' option filter alignments with low mapping quality.
 * Issue warnning and exit if no clipped reads found. 


## RSeQC v2.4
rewrite "geneBody_coverage.py"

* Memory-efficient: consumed < 100M ## RAM
* Flexible input to handle one or more BAM files::

 * Input a singel BAM file. eg: **-i test.bam**
 * Input several BAM files (separated by ","). eg: **-i test1.bam,test2.bam,test3.bam**
 * Input plain text file containing the path of BAM file(s). eg: **-i bam_path.txt**
 * Input a directory containing BAM file(s). eg: **-i my_folder**

* Generate heatmap to visualize gene body coverage over many samples.

## RSeQC v2.3.9
 
* Add bam2fq.py. Transform BAM files into fastq format.
* bugs fixed. 

## RSeQC v2.3.7

* bam_stat.py: Now counts 'Proper-paired reads map to different chrom'
* bam2wig.py: automatically call 'wigToBigwig' if it can be found in system $PATH
* inner_distance.py: add 'PE_within_diff_chrom'

## RSeQC v2.3.3

* Minor bugs fixed. 

## RSeQC v2.3.2

* Add split_bam.py: Split orignal BAM file into small BAM files based on provided gene list. User can use this module to estimate ribosome ## RNA amount if the input gene list is ribosomal ## RNA.
* Add  read_hexamer.py: Calculate hexamer frequency for multiple input files (fasta or fastq).
* Some parts were optimized and runs little faster.

## RSeQC v2.3.1

* Add normalization option to bam2wig.py. With this option, user can normalize different sequencing depth into the same scale when converting BAM into wiggle format.
* Add another script. geneBody_coverage2.py. This script uses BigWig? instead of BAM as input, and requires much less memory (~ 200M) 
