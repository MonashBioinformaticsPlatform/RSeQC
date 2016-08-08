#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
For each gene, check whether the RPKM value was saturated or not. Saturation analysis is based on 
re-sampling. For example, sample 5%, 10%, ... , 95%, 100% from total mapped reads, then 
calculate RPKM value for each step. Strand specific sequencing protocol is supported.
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
import numpy as np
from time import strftime
import operator
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

def normalize(lst):
	'''normalize all numbers between 0 and 1'''
	norm_lst=[]
	if max(lst) == min(lst):
		return norm_lst
	if max(lst) - min(lst)==0:
		return norm_lst
	for i in lst:
		norm_lst.append(  (i-min(lst))/(max(lst)-min(lst))  )  
	return norm_lst

def square_error(lst):
	'''transform list into normalized squared error (squared error divided by range)'''
	SE = []
	true_rpkm = lst[-1]
	rang = max(lst) - min(lst)
	if true_rpkm == 0:
		return None
	if rang == 0:
		return None
	for i in lst:
		SE.append(  abs(i - true_rpkm)/true_rpkm  )  
	return SE


def show_saturation (infile,outfile,rpkm_cut=0.01):
	
	RPKM_values = collections.defaultdict(list)
	RPKM_mean = {}
	gene_count = 0
	Quan = {'Q1':[0,0.25],'Q2':[0.25,0.5],'Q3':[0.5,0.75],'Q4':[0.75,1]}
	ROUT = open(outfile,'w')
	
	for line in open(infile):
		line=line.strip()
		fields=line.split()
		if fields[0].startswith('#'):
			head = [i.replace('%','') for i in fields[6:]]
			continue
		mykey = '\t'.join(fields[0:6])
		myvalue = [float(i) for i in fields[6:]]
		if max(myvalue) == 0: continue
		if max(myvalue) - min(myvalue) == 0: continue
		if np.mean(myvalue) < rpkm_cut: continue
		
		RPKM_values[mykey] = square_error(myvalue)
		RPKM_mean[mykey] = np.mean(myvalue)
		gene_count += 1
		if (len(head)==0):
			print >>sys.stderr, "No head line found, exit."
			sys.exit(1)
	
	
	print >>ROUT, "pdf('%s')" % (outfile.replace('.r','.pdf'))
	print >>ROUT, 'par(mfrow=c(2,2))'
	for quantile in sorted(Quan):
		line_count = 0
		norm_RPKM = collections.defaultdict(list)
		for k,v in sorted (RPKM_mean.iteritems(), key=operator.itemgetter(1)):	
			line_count += 1
			if (line_count > gene_count * Quan[quantile][0]) and (line_count <= gene_count * Quan[quantile][1]):
				for i,j in enumerate(RPKM_values[k]):
					norm_RPKM[head[i]].append(str(j))
		print >>ROUT,"name=c(%s)" % (','.join(head[:-1]))
		for i in head[:-1]:
			print >>ROUT, "S%s=c(%s)" % (i, ','.join(norm_RPKM[i]))
		print >>ROUT, "boxplot(%s,names=name,outline=F,ylab='Percent Relative Error',main='%s',xlab='Resampling percentage')" % (','.join(['100*S' + i for i in head[:-1]]),quantile)
	print >>ROUT, 'dev.off()'
				

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format. [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat. [required]")
	parser.add_option("-d","--strand",action="store",type="string",dest="strand_rule",default=None,help="How read(s) were stranded during sequencing. For example: --strand='1++,1--,2+-,2-+' means that this is a pair-end, strand-specific RNA-seq, and the strand rule is: read1 mapped to '+' => parental gene on '+'; read1 mapped to '-' => parental gene on '-'; read2 mapped to '+' => parental gene on '-'; read2 mapped to '-' => parental gene on '+'.  If you are not sure about the strand rule, run \'infer_experiment.py' default=%default (Not a strand specific RNA-seq data)")
	parser.add_option("-l","--percentile-floor",action="store",type="int",dest="percentile_low_bound",default=5, help="Sampling starts from this percentile. A integer between 0 and 100. default=%default")
	parser.add_option("-u","--percentile-ceiling",action="store",type="int",dest="percentile_up_bound",default=100, help="Sampling ends at this percentile. A integer between 0 and 100. default=%default")
	parser.add_option("-s","--percentile-step",action="store",type="int",dest="percentile_step",default=5, help="Sampling frequency. Smaller value means more sampling times. A integer between 0 and 100. default=%default")	
	parser.add_option("-c","--rpkm-cutoff",action="store",type="float",dest="rpkm_cutoff",default=0.01, help="Transcripts with RPKM smaller than this number will be ignored in visualization plot. default=%default")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")
	
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file):
		parser.print_help()
		sys.exit(0)
	if options.percentile_low_bound <0 or options.percentile_low_bound >100:
		print >>sys.stderr, "percentile_low_bound must be larger than 0 and samller than 100"
		sys.exit(0)
	if options.percentile_up_bound <0 or options.percentile_up_bound >100:
		print >>sys.stderr, "percentile_up_bound must be larger than 0 and samller than 100"
		sys.exit(0)
	if options.percentile_up_bound < options.percentile_low_bound:
		print >>sys.stderr, "percentile_up_bound must be larger than percentile_low_bound"
		sys.exit(0)
	if options.percentile_step <0 or options.percentile_step > options.percentile_up_bound:
		print >>sys.stderr, "percentile_step must be larger than 0 and samller than percentile_up_bound"
		sys.exit(0)
	if os.path.exists(options.input_file):
		obj = SAM.ParseBAM(options.input_file)
		obj.saturation_RPKM(outfile=options.output_prefix, refbed=options.refgene_bed, sample_start=options.percentile_low_bound,sample_end=options.percentile_up_bound,sample_step=options.percentile_step,strand_rule=options.strand_rule, q_cut  = options.map_qual)
		show_saturation(infile=options.output_prefix + ".eRPKM.xls", outfile=options.output_prefix + ".saturation.r",rpkm_cut = options.rpkm_cutoff)
		try:
			subprocess.call("Rscript " + options.output_prefix + ".saturation.r", shell=True)
		except:
			pass
	else:
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)
		


if __name__ == '__main__':
	main()
