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
- [Contact](#contact)
- [Reference](#reference)

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

## Installation
 
You will need either `sudo` or [virtualenvs](ttp://docs.python-guide.org/en/latest/dev/virtualenvs/) (which is my preferred method). If you are you going to use `sudo` please prefix `python setup.py install` and `pip install numpy` with `sudo`.

```BASH
git clone --branch fresh https://github.com/MonashBioinformaticsPlatform/RSeQC.git
cd RSeQC
python setup.py install
rseqc --help
```

I haven't figured why, but `numpy` needs to be installed separately. It doesn't get pulled correctly from the dependencies list in `setup.up`.

```BASH
pip install numpy
```

## Input format

- [BED](http://genome.ucsc.edu/FAQ/FAQformat.html) file is tab separated, 12-column, plain text file to represent gene models
- [GTF](http://mblab.wustl.edu/GTF22.html) file is also represents gene models. This is an alternative file to BED12
- [SAM/BAM](http://www.htslib.org/doc/sam.html) file holds information about read alignment to the reference genome. 

## Contact

- Liguo Wang: wangliguo78@gmail.com
- Shengqin Wang: wzsqwang@gmail.com
- Wei Li: superliwei@gmail.com 

## Reference

- Wang, L., Wang, S., & Li, W. (2012). **RSeQC: quality control of RNA-seq experiments**. *Bioinformatics* (Oxford, England), 28(16), 2184–2185. http://doi.org/10.1093/bioinformatics/bts356
- Wang, L., Nie, J., Sicotte, H., Li, Y., Eckel-Passow, J. E., Dasari, S., et al. (2016). **Measure transcript integrity using RNA-seq data**. *BMC Bioinformatics*, 17(1), 1–16. http://doi.org/10.1186/s12859-016-0922-z
