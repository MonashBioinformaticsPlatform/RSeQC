#!/usr/bin/env python
'''Manipulate two bigwig files'''
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import string
from optparse import OptionParser

from bx.bbi.bigwig_file import BigWigFile
from qcmodule import BED
from qcmodule import twoList
import numpy
__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"

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
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	
	parser.add_option("-i","--bwfile1",action="store",type="string",dest="BigWig_File1",help="One BigWig file")
	parser.add_option("-j","--bwfile2",action="store",type="string",dest="BigWig_File2",help="Another BigWig file")
	parser.add_option("-a","--action",action="store",type="string",dest="action",help='After pairwise align two bigwig files, perform the follow actions (Only select one keyword):"Add" = add signals. "Average" = average signals. "Division"= divide bigwig2 from bigwig1. Add 1 to both bigwig. "Max" = pick the signal that is larger. "Min" = pick the signal that is smaller. "Product" = multiply signals. "Subtract" = subtract signals in 2nd bigwig file from the corresponiding ones in the 1st bigwig file. "geometricMean" = take the geometric mean of signals.')
	parser.add_option("-o","--output",action="store",type="string",dest="output_wig",help="Output wig file")
	parser.add_option("-s","--chromSize",action="store",type="string",dest="chromSize",help="Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name, second column is size of the chromosome.")
	parser.add_option("-c","--chunk",action="store",type="int",dest="chunk_size",default=100000,help="Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=%default (bp)")
	parser.add_option("-m","--min_signal",action="store",type="float",dest="min_score",default=0.0,help="To redude the size of output wigfile, genomic positions with signal value smaller than (<) this threshold will be filtered out. default=%default")
	(options,args)=parser.parse_args()
	
	if not (options.BigWig_File1 and options.BigWig_File2  and options.output_wig and options.chromSize):
		parser.print_help()
		sys.exit(0)
	OUT=open(options.output_wig,'w')
	bw1 = BigWigFile( file=open(options.BigWig_File1) )
	bw2 = BigWigFile( file=open(options.BigWig_File2) )
	chrom_sizes = load_chromsize(options.chromSize)
	for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
		print >>sys.stderr, "Processing " + chr_name + " ..."
		OUT.write('variableStep chrom='+chr_name+'\n')
		for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
			coord = interval[1]
			try:
				bw_signal1 = bw1.get_as_array(chr_name,interval[1],interval[2])
			except:
				bw_signal1 = numpy.array()
			try:
				bw_signal2 = bw2.get_as_array(chr_name,interval[1],interval[2])
			except:
				bw_signal2 = numpy.array()
			if bw_signal1 is None and bw_signal2 is None:
				continue
			if numpy.isnan(numpy.nansum(bw_signal1)) and numpy.isnan(numpy.nansum(bw_signal2)):
				continue
			if len(bw_signal1) == 0 and len(bw_signal2) == 0:
				continue
			bw_signal1 = numpy.nan_to_num( bw_signal1 )
			bw_signal2 = numpy.nan_to_num( bw_signal2 )
		
			call_back = getattr(twoList,options.action)
			for v in call_back(bw_signal1,bw_signal2):
				coord +=1
				if v >= options.min_score: print >>OUT, "%d\t%.2f" % (coord,v)

if __name__=='__main__':
	main()
