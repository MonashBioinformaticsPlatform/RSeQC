#!/usr/bin/env python
'''
Calculate the RNA-seq reads coverage over gene body.
This module uses bigwig file as input.
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
import numpy as np
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.bbi.bigwig_file import BigWigFile

#import my own modules
from qcmodule import SAM
from qcmodule import mystat
#changes to the paths
__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def coverageGeneBody_bigwig(bigFile,refbed,outfile,gtype="png"):
	'''Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divided
	into 100 regsions. bigFile is bigwig format file'''
	if refbed is None:
		print >>sys.stderr,"You must specify a bed file representing gene model\n"
		exit(0)
	OUT1 = open(outfile + ".geneBodyCoverage_plot.r",'w')
	OUT2 = open(outfile + ".geneBodyCoverage.txt",'w')
	
	bw = BigWigFile( file = open(bigFile) )
	print >>sys.stderr, "calculating coverage over gene body ..."
	coverage=collections.defaultdict(int)
	flag=0
	gene_count = 0
	for line in open(refbed,'r'):
		try:
			if line.startswith(('#','track','browser')):continue  
			gene_count += 1
           	# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
				
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue
		gene_all_base=[]
		percentile_base=[]
		mRNA_len =0
		flag=0
		for st,end in zip(exon_starts,exon_ends):
			gene_all_base.extend(range(st+1,end+1))		#0-based coordinates on genome
			mRNA_len = len(gene_all_base)
			if mRNA_len <100:
				flag=1
				break
		if flag==1: continue
		if strand == '-':
			gene_all_base.sort(reverse=True)			#deal with gene on minus stand
		else:
			gene_all_base.sort(reverse=False)
		percentile_base = mystat.percentile_list (gene_all_base)	#get 101 points from each gene's coordinates
			
		for i in range(0,len(percentile_base)):
			#try:
			sig = bw.get_as_array(chrom,percentile_base[i]-1,percentile_base[i])
			if sig is None:continue
			coverage[i] += np.nan_to_num(sig[0])
			#except:
			#	continue
		print >>sys.stderr, "  %d genes finished\r" % gene_count,

	x_coord=[]
	y_coord=[]
	print >>OUT2, "percentile\tcount"
	for i in coverage:
		x_coord.append(str(i))
		y_coord.append(str(coverage[i]))
		print >>OUT2, str(i) + '\t' + str(coverage[i])
		
	print >>OUT1, "%s(\'%s\')" % (gtype, outfile + ".geneBodyCoverage." + gtype)
	print >>OUT1, "x=1:100"
	print >>OUT1, "y=c(" + ','.join(y_coord) + ')'
	print >>OUT1, "plot(x,y/%s,xlab=\"percentile of gene body (5'->3')\",ylab='average wigsum',type='s')" % gene_count
	print >>OUT1, "dev.off()"
	
def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Coverage signal file in bigwig format")
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-t","--graph-type",action="store",type="string",dest="graph_type",default="pdf",help="Graphic file type in \"pdf\", \"jpeg\", \"bmp\", \"bmp\", \"tiff\" or \"png\".default=%default [optional]")
	(options,args)=parser.parse_args()

	gt = options.graph_type.lower()
	if gt not in ("pdf","png",'bmp','jpeg','tiff'):
		print >>sys.stderr, "graphic file type must be 'pdf' or 'png'"
		parser.print_help()
		sys.exit(0)
	if not (options.output_prefix and options.input_file and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print >>sys.stderr, '\n\n' + options.ref_gene_model + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)
	if os.path.exists(options.input_file):
		coverageGeneBody_bigwig(options.input_file,options.ref_gene_model,options.output_prefix,gtype=options.graph_type)
		try:
			subprocess.call("Rscript " + options.output_prefix + '.geneBodyCoverage_plot.r',shell=True)
		except:
			print >>sys.stderr, "Cannot generate plot from " + options.output_prefix + '.geneBodyCoverage_plot.r'
			pass
	else:
		print >>sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)




if __name__ == '__main__':
        main()
 
