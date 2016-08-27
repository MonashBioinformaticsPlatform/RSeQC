#!/usr/bin/env python
'''
Calculate the inner distance (insert size)  of RNA-seq fragments. 

               RNA fragment
 _________________||_________________
|                                    |
|                                    |
||||||||||------------------||||||||||
  read_1      insert_size     read_2

fragment size = read_1 + insert_size + read_2
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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s)")
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene",help="Reference gene model in BED format.")
	parser.add_option("-k","--sample-size",action="store",type="int",dest="sampleSize",default=1000000,help="Number of read-pairs used to estimate inner distance. default=%default")
	parser.add_option("-l","--lower-bound",action="store",type="int",dest="lower_bound_size",default=-250,help="Lower bound of inner distance (bp). This option is used for ploting histograme. default=%default")
	parser.add_option("-u","--upper-bound",action="store",type="int",dest="upper_bound_size",default=250,help="Upper bound of inner distance (bp). This option is used for plotting histogram. default=%default")
	parser.add_option("-s","--step",action="store",type="int",dest="step_size",default=5,help="Step size (bp) of histograme. This option is used for plotting histogram. default=%default")	
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")

	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file and options.ref_gene):
		parser.print_help()
		sys.exit(0)
	for input_file in ([options.input_file,options.ref_gene]):
		if not os.path.exists(input_file):
			print >>sys.stderr, '\n\n' + input_file + " does NOT exists" + '\n'
			parser.print_help()
			sys.exit(0)
	if options.step_size <=0:
		print >>sys.stderr, "step size is a positive interger"
		sys.exit(0)
	obj = SAM.ParseBAM(options.input_file)
	obj.mRNA_inner_distance(outfile=options.output_prefix,low_bound=options.lower_bound_size,up_bound=options.upper_bound_size,step=options.step_size,refbed=options.ref_gene,sample_size=options.sampleSize, q_cut = options.map_qual)
	try:
		subprocess.call("Rscript " + options.output_prefix + '.inner_distance_plot.r',shell=True)
	except:
		print >>sys.stderr, "Cannot generate pdf file from " + options.output_prefix + '.inner_distance_plot.r'
		pass

if __name__ == '__main__':
	main()
