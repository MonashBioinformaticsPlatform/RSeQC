#!/usr/bin/env python
'''
Convert BAM file into wig file. BAM file must be sorted and indexed using SAMtools.
Note: SAM format file is not supported.
'''

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

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
from qcmodule import BED
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('class.log','a')
	print >>sys.stderr,mesg
	print >>LOG,mesg

def load_chromsize(file):
	'''read chrom.size file'''
	chromSize={}
	for line in open(file,'r'):
		if line.startswith('#'):continue
		if not line.strip():continue
		fields = line.strip().split()
		chromSize[fields[0]] = int(fields[1])
	return chromSize

			
def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM format. BAM file must be sorted and indexed using samTools. .bam and .bai files should be placed in the same directory.")
	parser.add_option("-s","--chromSize",action="store",type="string",dest="chromSize",help="Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name/ID, second column is chromosome size. Chromosome name (such as \"chr1\") should be consistent between this file and the BAM file.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output wiggle files(s). One wiggle file will be generated for non strand-specific data, two wiggle files (\"Prefix_Forward.wig\" and \"Prefix_Reverse.wig\") will be generated for strand-specific RNA-seq data.")
	parser.add_option("-t","--wigsum",action="store",type="int",dest="total_wigsum",help="Specified wigsum. Eg: 1,000,000,000 equals to coverage of 10 million 100nt reads. Ignore this option to disable normalization")
	parser.add_option("-u","--skip-multi-hits",action="store_true",dest="skip_multi",help="Skip non-unique hit reads.")
	parser.add_option("-d","--strand",action="store",type="string",dest="strand_rule",default=None,help="How read(s) were stranded during sequencing. For example: --strand='1++,1--,2+-,2-+' means that this is a pair-end, strand-specific RNA-seq data, and the strand rule is: read1 mapped to '+' => parental gene on '+'; read1 mapped to '-' => parental gene on '-'; read2 mapped to '+' => parental gene on '-'; read2 mapped to '-' => parental gene on '+'.  If you are not sure about the strand rule, run \'infer_experiment.py' default=%default (Not a strand specific RNA-seq data).")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality to determine \"uniquely mapped\". default=%default")

	(options,args)=parser.parse_args()
	
	
	if not (options.output_prefix and options.input_file and options.chromSize and options.output_prefix):
		parser.print_help()
		sys.exit(0)
	for file in (options.input_file,options.chromSize):
		if not os.path.exists(file):
			print >>sys.stderr, '\n\n' + file + " does NOT exists" + '\n'
			sys.exit(0)
	if not os.path.exists(options.input_file + '.bai'):
		print >>sys.stderr, "index file " + options.input_file + '.bai' + " does not exists"
		sys.exit(0)
	
	if options.skip_multi:print "Skip multi-hits:True"
	else:print "Skip multi-hits:False"


	chromSizes = load_chromsize(options.chromSize)
	
	norm_factor=None
	if options.total_wigsum:
		obj = SAM.ParseBAM(options.input_file)
		wig_sum = obj.calWigSum(chrom_sizes = chromSizes, skip_multi=options.skip_multi)
		print >>sys.stderr, "\n\ntotal wigsum is:" + str(wig_sum) + '\n'
		try:
			norm_factor = options.total_wigsum / wig_sum
		except:
			norm_factor = None
			
	obj = SAM.ParseBAM(options.input_file)		
	obj.bamTowig(outfile = options.output_prefix, chrom_sizes = chromSizes, chrom_file = options.chromSize, q_cut = options.map_qual, skip_multi=options.skip_multi,strand_rule = options.strand_rule, WigSumFactor=norm_factor)

if __name__ == '__main__':
	main()
