# RSeQC: An RNA-seq Quality Control Package

> This is a fork of original RSeQC package from [sorceforge site](http://rseqc.sourceforge.net)
> I'm making drastic rearrangement to this package to make it easier to follow. I'm also making changes to the code base
> At this stage only `read_dist` (read_distribution) and `bam_stats` (bam_stat) modules have been incorporated and both now can be accessed from main executable `scripts/rseqc` 

**Original message**

> RSeQC package provides a number of useful modules that can comprehensively evaluate high
> throughput sequence data especially RNA-seq data. Some basic modules quickly inspect sequence
> quality, nucleotide composition bias, PCR bias and GC bias, while RNA-seq specific modules
> evaluate sequencing saturation, mapped reads distribution, coverage uniformity, strand specificity, transcript level RNA integrity etc.

## Table of content

- [Quick start](#quick-start)
- [Installation](#installation)
- [Input formats](#input-format)

## Quick start

Once installed use main executable file `rseqc` to run any of the sub-commands (modules)

e.g

```BASH
rseqc read_dist --input_file yourBamFile.bam --gene_models yourGTFfile.gtf
```

OR

```BASH
rseqc read_dist --input_file yourBamFile.bam --gene_models yourBED12file.bed --file_type bed
```


Install RSeQC from source code (Not recommended)
----------------------------------------------------
 
Prerequisite: gcc; `python2.7 <http://www.python.org/getit/releases/2.7/>`_;  `numpy <http://numpy.scipy.org/>`_; `R <http://www.r-project.org/>`_
 
Install RSeQC (Example)::
 
 tar zxf RSeQC-VERSION.tar.gz
  
 cd RSeQC-VERSION
 
 #type 'python setup.py install --help' to see options
 python setup.py install	#Note this requires root privilege
 or
 python setup.py install --root=/home/user/XXX/		#install RSeQC to user specificed location, does NOT require root privilege
 
 #This is only an example. Change path according to your system configuration
 export PYTHONPATH=/home/user/lib/python2.7/site-packages:$PYTHONPATH 
 
 #This is only an example. Change path according to your system configuration
 export PATH=/home/user/bin:$PATH

Finally, type: python -c 'from qcmodule import SAM'. If no error message comes out, RSeQC
modules have been installed successfully. 

## Input format

RSeQC accepts 4 file formats as input:

- `BED <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ file is tab separated, 12-column, plain text file to represent gene model. Here is an `example <http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/sample.bed>`_.
- `SAM <http://samtools.sourceforge.net/>`_ or `BAM <http://genome.ucsc.edu/goldenPath/help/bam.html>`_ files are used to store reads alignments. SAM file is human readable plain text file, while BAM is  binary version of SAM, a compact and index-able representation of reads alignments. Here is an `example <http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/sample.sam>`_.
- Chromosome size file is a two-column, plain text file. Here is an `example <http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/sample.hg19.chrom.sizes>`_ for human hg19 assembly. Use this `script <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes>`_ to download chromosome size files of other genomes.
- `Fasta <http://en.wikipedia.org/wiki/FASTA_format>`_ file.

## Fetch chromsome size file from UCSC

download this `script <http://sourceforge.net/projects/rseqc/files/other/fetchChromSizes/download>`_
and save as 'fetchChromSizes'::
 
 # Make sure it's executable 
 chmod +x fetchChromSizes
 
 fetchChromSizes hg19 >hg19.chrom.sizes
 
 fetchChromSizes danRer7  >zebrafish.chrom.sizes

## Contact

- Liguo Wang: wangliguo78@gmail.com
- Shengqin Wang: wzsqwang@gmail.com
- Wei Li: superliwei@gmail.com 

## Reference

- Wang, L., Wang, S., & Li, W. (2012). **RSeQC: quality control of RNA-seq experiments**. *Bioinformatics* (Oxford, England), 28(16), 2184–2185. http://doi.org/10.1093/bioinformatics/bts356
- Wang, L., Nie, J., Sicotte, H., Li, Y., Eckel-Passow, J. E., Dasari, S., et al. (2016). **Measure transcript integrity using RNA-seq data**. *BMC Bioinformatics*, 17(1), 1–16. http://doi.org/10.1186/s12859-016-0922-z
