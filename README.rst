.. RSeQC documentation master file, created by
   sphinx-quickstart on Thu Aug  9 21:50:37 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   
RSeQC: An RNA-seq Quality Control Package
===========================================

RSeQC package provides a number of useful modules that can comprehensively evaluate high
throughput sequence data especially RNA-seq data. Some basic modules quickly inspect sequence
quality, nucleotide composition bias, PCR bias and GC bias, while RNA-seq specific modules
evaluate sequencing saturation, mapped reads distribution, coverage uniformity, strand specificity, transcript level RNA integrity etc.


Download
===============================

Download RSeQC
--------------------------------
* `RSeQC (v2.6.4) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.4.tar.gz/download>`_ (Note: Downloading "RSeQC-2.6.4.tar.gz" to local computer is unnecessary if you use **pip install RSeQC**)
* `RSeQC (v2.6.3) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.3.tar.gz/download>`_
* `RSeQC (v2.6.2) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.2.tar.gz/download>`_
* `RSeQC (v2.6.1) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.1.tar.gz/download>`_
* `RSeQC (v2.6) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.tar.gz/download>`_
* `RSeQC (v2.5) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.5.tar.gz/download>`_
* `RSeQC (v2.4) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.4.tar.gz/download>`_
* `RSeQC (v2.3.9) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.3.9.tar.gz/download>`_
* `RSeQC (v2.3.8) <http://sourceforge.net/projects/rseqc/files/RSeQC-2.3.8.tar.gz/download>`_
* `RSeQC (v2.3.7) <https://sourceforge.net/projects/rseqc/files/RSeQC-2.3.7.tar.gz/download>`_
* `RSeQC (v2.3.6) <http://rseqc.googlecode.com/files/RSeQC-2.3.6.tar.gz>`_
* `RSeQC (v2.3.3) <http://rseqc.googlecode.com/files/RSeQC-2.3.3.tar.gz>`_
* `RSeQC (v2.3.2) <http://rseqc.googlecode.com/files/RSeQC-2.3.2.tar.gz>`_
* `RSeQC (v2.3.1) <http://rseqc.googlecode.com/files/RSeQC-2.3.1.tar.gz>`_

Download test datasets
--------------------------------
 
Pair-end strand specific (Illumina). BAM file md5sum=fbd1fb1c153e3d074524ec70e6e21fb9

* `Pairend_StrandSpecific_51mer_Human_hg19.bam <https://drive.google.com/file/d/0BwAUopGWU_khNmozSHhWdDVncXc/view?usp=sharing>`_
* `Pairend_StrandSpecific_51mer_Human_hg19.bam.bai <https://drive.google.com/file/d/0BwAUopGWU_khc2g4akJlN25KVzQ/view?usp=sharing>`_
 
Pair-end  non-strand specific (Illumina). BAM file md5sum=ba014f6b397b8a29c456b744237a12de

* `Pairend_nonStrandSpecific_36mer_Human_hg19.bam <https://drive.google.com/file/d/0BwAUopGWU_khbjJPX3BxRzNtOWs/view?usp=sharing>`_
* `Pairend_nonStrandSpecific_36mer_Human_hg19.bam.bai <https://drive.google.com/file/d/0BwAUopGWU_khUi02dGc0VjhORlk/view?usp=sharing>`_
  
Single-end strand specific (SOLiD). BAM file md5sum=b39951a6ba4639ca51983c2f0bf5dfce

* `SingleEnd_StrandSpecific_50mer_Human_hg19.bam <https://drive.google.com/file/d/0BwAUopGWU_khUDNTRHhFc29RQms/view?usp=sharing>`_
* `SingleEnd_StrandSpecific_50mer_Human_hg19.bam.bai <https://drive.google.com/file/d/0BwAUopGWU_khSVFLdUhGRjJpS0k/view?usp=sharing>`_
 
Download gene models (update on 08/07/2014)
--------------------------------------------


* `Human (hg38, hg19) <https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/>`_
* `Mouse (mm10,mm9) <https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/>`_
* `Fly (dm3) <https://sourceforge.net/projects/rseqc/files/BED/fly_D.melanogaster/>`_
* `Zebrafish (danRer7) <https://sourceforge.net/projects/rseqc/files/BED/Zebrafish_Danio_rerio/>`_

**NOTE:**

* BED file for other species and the most recent release of these files can be downloaded from `UCSC Table Browser <https://genome.ucsc.edu/cgi-bin/hgTables?command=start>`_ 
* Make sure the annotation file and the genome assembly are matched. For example, if you aligned RNA-seq reads to `hg19/GRCh37 <http://www.ncbi.nlm.nih.gov/assembly/2758/>`_ you should download `hg19/GRCh37 <http://www.ncbi.nlm.nih.gov/assembly/2758/>`_ based BED files. 


Download ribosome RNA (update on 07/08/2015)
---------------------------------------------
We only provide rRNA bed files for human and mouse. We download these ribosome RNAs from UCSC table browser,
we provide them here to facilitate users with NO WARRANTY in completeness.

* `GRCh38_rRNA.bed <http://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/GRCh38_rRNA.bed.gz/download>`_
* `hg19_rRNA.bed <http://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_rRNA.bed.gz/download>`_
* `mm10_rRNA.bed <http://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_rRNA.bed.gz/download>`_
* `mm9_rRNA.bed <http://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm9_rRNA.bed.gz/download>`_


Installation
===================================

Use pip to install RSeQC 
------------------------------------------

::

 pip install RSeQC
 
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

Input format
=============================

RSeQC accepts 4 file formats as input:

* `BED <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ file is tab separated, 12-column, plain text file to represent gene model. Here is an `example <http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/sample.bed>`_.
* `SAM <http://samtools.sourceforge.net/>`_ or `BAM <http://genome.ucsc.edu/goldenPath/help/bam.html>`_ files are used to store reads alignments. SAM file is human readable plain text file, while BAM is  binary version of SAM, a compact and index-able representation of reads alignments. Here is an `example <http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/sample.sam>`_.
* Chromosome size file is a two-column, plain text file. Here is an `example <http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/sample.hg19.chrom.sizes>`_ for human hg19 assembly. Use this `script <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes>`_ to download chromosome size files of other genomes.
* `Fasta <http://en.wikipedia.org/wiki/FASTA_format>`_ file.

**NOTE:**
If you have GFF/GTF format gene files, we found this `Perl script <https://code.google.com/p/ea-utils/source/browse/trunk/clipper/gtf2bed>`_ might be useful to convert them to BED. 

Fetch chromsome size file from UCSC
======================================
download this `script <http://sourceforge.net/projects/rseqc/files/other/fetchChromSizes/download>`_
and save as 'fetchChromSizes'::
 
 # Make sure it's executable 
 chmod +x fetchChromSizes
 
 fetchChromSizes hg19 >hg19.chrom.sizes
 
 fetchChromSizes danRer7  >zebrafish.chrom.sizes

Contact
================================
* Liguo Wang: wangliguo78@gmail.com
* Shengqin Wang: wzsqwang@gmail.com
* Wei Li: superliwei@gmail.com 

Reference
================================
* Wang, L., Wang, S., & Li, W. (2012). **RSeQC: quality control of RNA-seq experiments**. *Bioinformatics* (Oxford, England), 28(16), 2184–2185. http://doi.org/10.1093/bioinformatics/bts356
* Wang, L., Nie, J., Sicotte, H., Li, Y., Eckel-Passow, J. E., Dasari, S., et al. (2016). **Measure transcript integrity using RNA-seq data**. *BMC Bioinformatics*, 17(1), 1–16. http://doi.org/10.1186/s12859-016-0922-z
