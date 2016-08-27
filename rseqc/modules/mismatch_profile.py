#!/usr/bin/env python
'''
Calculate the distribution of mismatches across reads. Note that the "MD" tag must exist in BAM file.
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import string
import subprocess
from optparse import OptionParser
from time import strftime


#import my own modules
from qcmodule import SAM
from qcmodule import fasta
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
	print >>sys.stderr,mesg


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_bam",help='Input BAM file. [required]')
	parser.add_option("-l","--read-align-length",action="store",type="int", dest="read_alignment_length",help="Alignment length of read. It is usually set to the orignial read length. For example, all these cigar strings (\"101M\", \"68M140N33M\", \"53M1D48M\") suggest the read alignment length is 101. [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-n","--read-num",action="store",type="int",default=1000000, dest="read_number",help="Number of aligned reads with mismatches used to calculate the mismatch profile. default=%default")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality. default=%default")
	(options,args)=parser.parse_args()

	if not (options.input_bam):
		parser.print_help()
		sys.exit(0)
	for f in ([options.input_bam]):
		if not os.path.exists(f):
			print >>sys.stderr, '\n\n' + f + " does NOT exists" + '\n'
			parser.print_help()
			sys.exit(0)

	if not (options.output_prefix):
		print >>sys.stderr, '\n\n You must specify the output prefix'
		parser.print_help()
		sys.exit(0)
	
	if not (options.read_alignment_length):
		print >>sys.stderr, '\n\n You must specify read alignment length. It is usually the read length.'
		parser.print_help()
		sys.exit(0)	
		
	obj = SAM.ParseBAM(options.input_bam)
	obj.mismatchProfile(read_length = options.read_alignment_length, read_num = options.read_number, q_cut = options.map_qual,outfile = options.output_prefix)
	try:
		subprocess.call("Rscript " + options.output_prefix + '.mismatch_profile.r',shell=True)
	except:
		print >>sys.stderr, "Cannot generate pdf file from " + options.output_prefix + '.mismatch_profile.r'
		pass
	
if __name__ == '__main__':
	main()
