#!/usr/bin/env python
'''
Description: Convert alignments in BAM or SAM format into fastq format.
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import subprocess

from time import strftime

#import my own modules
from qcmodule import SAM
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



def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output fastq files(s).")
	parser.add_option("-s","--single-end",action="store_true",dest="single", help="Specificy '-s' or '--single-end' for single-end sequencing.")
	parser.add_option("-c","--compress",action="store_true",dest="gzip", help="Specificy '-c' or '--compress' to compress output fastq file(s) using 'gzip' command.")
	(options,args)=parser.parse_args()
	
	
	#print options.single
	if not (options.output_prefix and options.input_file):
		parser.print_help()
		sys.exit(0)
	if os.path.exists(options.input_file):
		obj = SAM.ParseBAM(options.input_file)
		if options.single is True:
			obj.bam2fq(prefix=options.output_prefix, paired = False)
			if options.gzip is True:
				try:
					print >>sys.stderr, "run gzip ... ",
					subprocess.call("gzip " + options.output_prefix + '.fastq', shell=True)
					print >>sys.stderr, "Done."
				except:
					pass
		else:
			obj.bam2fq(prefix=options.output_prefix, paired = True)
			if options.gzip is True:
				try:
					print >>sys.stderr, "run gzip ..."
					subprocess.call("gzip " + options.output_prefix + '.R1.fastq', shell=True)
					subprocess.call("gzip " + options.output_prefix + '.R2.fastq', shell=True)
					print >>sys.stderr, "Done."
				except:
					pass
	else:
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)

if __name__ == '__main__':
	main()
