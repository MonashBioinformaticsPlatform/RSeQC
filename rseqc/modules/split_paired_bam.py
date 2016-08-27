#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Split bam file (pair-end) into 2 single-end bam file
-------------------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import string
from optparse import OptionParser
import warnings
import string
import collections
import sets

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
import pysam

#import my own modules
from qcmodule import BED
from qcmodule import SAM
from qcmodule import bam_cigar

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"

	
def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output BAM files. \"prefix.R1.bam\" file contains the 1st read, \"prefix.R2.bam\" file contains the 2nd read")
	(options,args)=parser.parse_args()
		
	if not (options.input_file):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file):
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		sys.exit(0)		
		
	samfile = pysam.Samfile(options.input_file,'rb')
	OUT1 = pysam.Samfile(options.output_prefix + '.R1.bam','wb',template=samfile)	#bam file containing reads hit to exon region
	OUT2 = pysam.Samfile(options.output_prefix + '.R2.bam','wb',template=samfile)	#bam file containing reads not hit to exon region
	OUT3 = pysam.Samfile(options.output_prefix + '.unmap.bam','wb',template=samfile)	#bam file containing reads not hit to exon region
	
	total_alignment = 0
	r1_alignment = 0
	r2_alignment = 0
	unmapped = 0
	
	print >>sys.stderr, "spliting " + options.input_file + " ...",
	try:
		while(1):
			new_alignment = pysam.AlignedRead()     # create AlignedRead object
			old_alignment = samfile.next()
			total_alignment += 1
			
			new_alignment.qname = old_alignment.qname      # 1st column. read name.
			#new_alignment.flag = old_alignment.flag        # 2nd column. subject to change. flag value 
			new_alignment.tid = old_alignment.tid          # 3rd column. samfile.getrname(tid) == chrom name
			new_alignment.pos = old_alignment.pos          # 4th column. reference Start position of the aligned part (of read) [0-based]
			new_alignment.mapq = old_alignment.mapq        # 5th column. mapping quality
			new_alignment.cigar= old_alignment.cigar       # 6th column. subject to change. 
			#new_alignment.rnext = old_alignment.rnext      # 7th column. tid of the reference (mate read mapped to)
			#new_alignment.pnext = old_alignment.pnext      # 8th column. position of the reference (0 based, mate read mapped to)
			#new_alignment.tlen = old_alignment.tlen        # 9th column. insert size
			new_alignment.seq = old_alignment.seq          # 10th column. read sequence. all bases.
			new_alignment.qual = old_alignment.qual        # 11th column. read sequence quality. all bases.
			new_alignment.tags = old_alignment.tags        # 12 - columns
			new_alignment.flag = 0x0000
			if old_alignment.is_unmapped:
				OUT3.write(old_alignment)
				unmapped += 1
				continue
			if old_alignment.is_reverse:
				new_alignment.flag = new_alignment.flag | 0x0010

			if old_alignment.is_secondary:
				new_alignment.flag = new_alignment.flag | 0x0100
			if old_alignment.is_qcfail:
				new_alignment.flag = new_alignment.flag | 0x0200
			if old_alignment.is_duplicate:
				new_alignment.flag = new_alignment.flag | 0x0400
			if old_alignment.is_read1:
				OUT1.write(new_alignment)
				r1_alignment += 1
			else:
				OUT2.write(new_alignment)
				r2_alignment += 1

	except StopIteration:
		print >>sys.stderr, "Done"
				
	print "%-55s%d" % ("Total records:",total_alignment)
	print "%-55s%d" % (options.output_prefix + 'Read 1:',r1_alignment)
	print "%-55s%d" % (options.output_prefix + 'Read 2:',r2_alignment)
	print "%-55s%d" % (options.output_prefix + 'Unmapped:',unmapped)
	
if __name__ == '__main__':
	main()
