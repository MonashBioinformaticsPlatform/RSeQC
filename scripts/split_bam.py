#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
Split bam file according to input gene list (bed).
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

def searchit(exon_range, exon_list):
	'''return 1 if find, return 0 if cannot find'''
	for chrom, st, end in exon_list:
		if chrom.upper() not in exon_range:
			return 0
		elif len(exon_range[chrom].find(st,end)) >=1:
			return 1
	return 0

def build_bitsets(list):
	'''build intevalTree from list'''
	ranges={}
	for l in list:
		chrom =l[0].upper()
		st = int(l[1])
		end = int(l[2])
		if chrom not in ranges:
			ranges[chrom] = Intersecter()
		ranges[chrom].add_interval( Interval( st, end ) )
	return ranges

	
def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed.")
	parser.add_option("-r","--genelist",action="store",type="string",dest="gene_list",help="Gene list in bed foramt. All reads hits to exon regions (defined by this gene list) will be saved into one BAM file, the remaining reads will saved into another BAM file.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output BAM files. \"prefix.in.bam\" file contains reads mapped to the gene list specified by \"-r\", \"prefix.ex.bam\" contains reads that cannot mapped to gene list. \"prefix.junk.bam\" contains qcfailed or unmapped reads.")
	(options,args)=parser.parse_args()
		
	if not (options.input_file and options.gene_list):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.gene_list):
		print >>sys.stderr, '\n\n' + options.gene_list + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file):
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		sys.exit(0)		
	
	#build bitset for gene list
	print >>sys.stderr, 'reading ' + options.gene_list + ' ... ',
	obj = BED.ParseBED(options.gene_list)
	exons = obj.getExon()
	exon_ranges = build_bitsets(exons)
	print >>sys.stderr, 'Done'
	
	samfile = pysam.Samfile(options.input_file,'rb')
	out1 = pysam.Samfile(options.output_prefix + '.in.bam','wb',template=samfile)	#bam file containing reads hit to exon region
	out2 = pysam.Samfile(options.output_prefix + '.ex.bam','wb',template=samfile)	#bam file containing reads not hit to exon region
	out3 = pysam.Samfile(options.output_prefix + '.junk.bam','wb',template=samfile)	#bam file containing reads not hit to exon region
	
	total_alignment = 0
	in_alignment = 0
	ex_alignment = 0
	bad_alignment = 0
	print >>sys.stderr, "spliting " + options.input_file + " ...",
	try:
		while(1):
			aligned_read = samfile.next()
			total_alignment += 1
			
			if aligned_read.is_qcfail:
				bad_alignment +=1
				out3.write(aligned_read)
				continue
			if aligned_read.is_unmapped:
				bad_alignment +=1
				out3.write(aligned_read)
				continue
			
			chrom = samfile.getrname(aligned_read.tid)
			chrom=chrom.upper()	
			read_start = aligned_read.pos
			mate_start = aligned_read.mpos
				
			#read_exons = bam_cigar.fetch_exon(chrom, aligned_read.pos, aligned_read.cigar)
			if aligned_read.mate_is_unmapped:	#only one end mapped
				if chrom not in exon_ranges:
					out2.write(aligned_read)
					ex_alignment += 1
					continue		
				else:		
					if len(exon_ranges[chrom].find(read_start, read_start +1)) >= 1:
						out1.write(aligned_read)
						in_alignment += 1
						continue
					elif len(exon_ranges[chrom].find(read_start, read_start +1)) == 0:
						out2.write(aligned_read)
						ex_alignment += 1
						continue
			else:							#both end mapped
				if chrom not in exon_ranges:
					out2.write(aligned_read)
					ex_alignment += 1
					continue
				else:
					if (len(exon_ranges[chrom].find(read_start, read_start +1)) >= 1) or (len(exon_ranges[chrom].find(mate_start, mate_start +1)) >= 1):
						out1.write(aligned_read)
						in_alignment += 1
					else:
						out2.write(aligned_read)
						ex_alignment += 1
				
	except StopIteration:
		print >>sys.stderr, "Done"
				
	print "%-55s%d" % ("Total records:",total_alignment)
	print "%-55s%d" % (options.output_prefix + '.in.bam (Reads consumed by input gene list):',in_alignment)
	print "%-55s%d" % (options.output_prefix + '.ex.bam (Reads not consumed by input gene list):',ex_alignment)
	print "%-55s%d" % (options.output_prefix + '.junk.bam (qcfailed, unmapped reads):',bad_alignment)
	
if __name__ == '__main__':
	main()
