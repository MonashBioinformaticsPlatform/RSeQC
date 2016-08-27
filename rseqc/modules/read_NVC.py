#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
For each nucleotide  position of read (5'->3'), check the nucleotide frequency. The generated R script will
gives NVC (Nucleotide Versus Cycle) plot.
-------------------------------------------------------------------------------------------------'''

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
import math
import sets
from time import strftime
import subprocess

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

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


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('class.log','a')
	print >>sys.stderr,mesg
	print >>LOG,mesg


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Input file in BAM or SAM format.[required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-x","--nx",action="store_true",dest="unknown_nucleotide",help="Flag option. Presense of this flag tells program to include N,X in output NVC plot [required]")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file):
		parser.print_help()
		sys.exit(0)
	if os.path.exists(options.input_file):
		obj = SAM.ParseBAM(options.input_file)
		obj.readsNVC(outfile=options.output_prefix,nx=options.unknown_nucleotide, q_cut = options.map_qual)
		try:
			subprocess.call("Rscript " + options.output_prefix +  ".NVC_plot.r",shell=True)
		except:
			pass
	else:
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)

if __name__ == '__main__':
	main()
