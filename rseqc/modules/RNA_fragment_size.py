#!/usr/bin/env python
'''
calculate fragment size for each gene/transcript. For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
'''

import pysam
import sys,os
from numpy import mean,median,std
from optparse import OptionParser

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"

"""
def overlap_length1(lst1, lst2):
	l = 0
	for chr, st, end  in BED.intersectBed3(lst1, lst2):
		l += end -st
	return l
"""
def overlap_length2(lst1, lst2):
	l = 0
	for x in lst1:
		for y in lst2:
			l += len(range(max(x[0], y[0]), min(x[-1], y[-1])+1))
	return l
		
	
def fragment_size(bedfile, samfile,qcut=30,ncut=5):
	'''calculate the fragment size for each gene'''
	for line in open(bedfile,'r'):
		exon_range = []	
		if line.startswith(('#','track','browser')):continue
		fields = line.split()
		chrom = fields[0]
		tx_start = int( fields[1] )
		tx_end = int( fields[2] )
		geneName = fields[3]
		trand    = fields[5].replace(" ","_")
		exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
		exon_starts = map((lambda x: x + tx_start ), exon_starts)
		exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
		exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
		geneID = "\t".join([str(i) for i in (chrom, tx_start, tx_end, geneName)])
		
		for st,end in zip(exon_starts,exon_ends):
			exon_range.append([st+1,end+1])
			#exon_range.append([chrom, st,end])
		
		try:
			alignedReads = samfile.fetch(chrom,tx_start,tx_end)
		except:
			yield '\t'.join([str(i) for i in (geneID, 0,0,0)])
			continue
		
		frag_sizes = []
		for aligned_read in alignedReads:
			if not aligned_read.is_paired:				#skip single sequencing
				continue
			if aligned_read.is_read2:
				continue
			if aligned_read.mate_is_unmapped:
				continue
			if aligned_read.is_qcfail:continue               #skip low quanlity                                      
			if aligned_read.is_duplicate:continue           #skip duplicate read
			if aligned_read.is_secondary:continue           #skip non primary hit           
			if aligned_read.mapq < qcut:
					continue
			
			read_st = aligned_read.pos
			mate_st = aligned_read.pnext
			if read_st > mate_st:
				(read_st, mate_st) = (mate_st, read_st)
			if read_st < tx_start or mate_st > tx_end:
				continue
			read_len = aligned_read.qlen
			map_range = [[read_st+1, mate_st]]
			#map_range = [[chrom, read_st, mate_st]]
			frag_len = overlap_length2(exon_range, map_range) + read_len
			frag_sizes.append(frag_len)
		if len(frag_sizes) < ncut:
			yield '\t'.join([str(i) for i in (geneID, len(frag_sizes), 0,0,0)])
		else:
			yield '\t'.join([str(i) for i in (geneID, len(frag_sizes), mean(frag_sizes), median(frag_sizes), std(frag_sizes))])


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_file",help="Input BAM file") 
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")
	parser.add_option("-n","--frag-num",action="store",type="int",dest="fragment_num",default=3,help="Minimum number of fragment. default=%default")

	(options,args)=parser.parse_args()

	if not (options.input_file and options.refgene_bed):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file + '.bai'):
		print >>sys.stderr, "cannot find index file of input BAM file"
		print >>sys.stderr, options.input_file + '.bai' + " does not exists"
		sys.exit(0)	

	for file in (options.input_file, options.refgene_bed):
		if not os.path.exists(file):
			print >>sys.stderr, file + " does NOT exists" + '\n'
			sys.exit(0)
	
	print '\t'.join([str(i) for i in ("chrom","tx_start", "tx_end","symbol","frag_count","frag_mean","frag_median","frag_std")])
	for tmp in fragment_size(options.refgene_bed, pysam.Samfile(options.input_file), options.map_qual, options.fragment_num):
		print tmp

if __name__ == '__main__':
        main()












			
			
