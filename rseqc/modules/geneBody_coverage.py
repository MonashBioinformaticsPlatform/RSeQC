#!/usr/bin/env python
'''
Calculate the RNA-seq reads coverage over gene body. 

Note:
1) Only input sorted and indexed BAM file(s). SAM format is not supported.
2) Genes/transcripts with mRNA length < 100 will be skipped (Number specified to "-l" cannot be < 100). 
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
import collections
import math
import sets
from time import strftime
import subprocess
from os.path import basename
import operator

#import third-party modules
from numpy import std,mean
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
import pysam

#import my own modules
from qcmodule import getBamFiles
from qcmodule import mystat
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


def valid_name(s):
	'''make sure the string 's' is valid name for R variable'''
	symbols = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'
	digit = '0123456789'
	rid = '_'.join(i for i in s.split())	#replace space(s) with '_'
	if rid[0] in digit:rid = 'V' + rid
	tmp = ''
	for i in rid:
		if i in symbols:
			tmp = tmp + i
		else:
			tmp = tmp + '_'
	return tmp
	

def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('log.txt','a')
	print >>sys.stderr,mesg
	print >>LOG,mesg


def pearson_moment_coefficient(lst):
	'''measure skewness'''
	mid_value = lst[int(len(lst)/2)]
	sigma = std(lst, ddof=1)
	tmp = []
	for i in lst:
		tmp.append(((i - mid_value)/sigma)**3)
	return mean(tmp)

def genebody_percentile(refbed, mRNA_len_cut = 100):
	'''
	return percentile points of gene body
	mRNA length < mRNA_len_cut will be skipped
	'''
	if refbed is None:
		print >>sys.stderr,"You must specify a bed file representing gene model\n"
		exit(0)
	
	g_percentiles = {}
	transcript_count = 0
	for line in open(refbed,'r'):
		try:
			if line.startswith(('#','track','browser')):continue  
			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
			geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
				
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
			transcript_count += 1
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue
		gene_all_base=[]
		mRNA_len =0
		flag=0
		for st,end in zip(exon_starts,exon_ends):
			gene_all_base.extend(range(st+1,end+1))		#1-based coordinates on genome
		if len(gene_all_base) < mRNA_len_cut:
			continue
		g_percentiles[geneID] = (chrom, strand, mystat.percentile_list (gene_all_base))	#get 100 points from each gene's coordinates
	printlog("Total " + str(transcript_count) + ' transcripts loaded')
	return g_percentiles

def genebody_coverage(bam, position_list):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	
	gene_finished = 0
	for chrom, strand, positions in position_list.values():
		coverage = {}
		for i in positions:
			coverage[i] = 0.0
		chrom_start = positions[0]-1
		if chrom_start <0: chrom_start=0
		chrom_end = positions[-1]
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
			
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions:
				continue
			if pileupcolumn.n == 0:
				coverage[ref_pos] = 0
				continue				
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read +=1
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
		if strand == '-':
			tmp = tmp[::-1]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]
		gene_finished += 1
		
		if gene_finished % 100 == 0:
			print >>sys.stderr, "\t%d transcripts finished\r" % (gene_finished),
	return 	aggreagated_cvg
		
def Rcode_write(dataset,file_prefix, format='pdf', colNum=100):
	'''generate R script for visualization'''
	ROUT = open(file_prefix + '.r','w')
	names=[]
	datas=[]
	for name, data, tmp in dataset:
		names.append(name)
		datas.append(data)
		print >>ROUT, name + ' <- c(' + ','.join([str(i) for i in data]) + ')'
		
	tick_pos = [1,10,20,30,40,50,60,70,80,90,100]
	tick_lab = [1,10,20,30,40,50,60,70,80,90,100]
	
	# do not generate heatmap if only 1 sample
	if len(names) >=3:
		print >>ROUT, 'data_matrix' + ' <- matrix(c(' + ','.join(names) + '), byrow=T, ' +  'ncol=' + str(colNum) + ')'
		print >>ROUT, 'rowLabel <- c(' + ','.join(['"' + i + '"' for i in names]) + ')'
		print >>ROUT, '\n'
		print >>ROUT, '%s(\"%s.%s\")' % (format.lower(),file_prefix + ".heatMap",format.lower())
		print >>ROUT, 'rc <- cm.colors(ncol(data_matrix))'
		print >>ROUT, 'heatmap(data_matrix' + ', scale=c(\"none\"),keep.dendro=F, labRow = rowLabel ' + ',Colv = NA,Rowv = NA,labCol=NA,col=cm.colors(256),margins = c(6, 8),ColSideColors = rc,cexRow=1,cexCol=1,xlab="Gene body percentile (5\'->3\')", add.expr=x_axis_expr <- axis(side=1,at=c(%s),labels=c(%s)))' % (','.join([str(i) for i in tick_pos]), ','.join(['"' + str(i) + '"' for i in tick_lab]))
		print >>ROUT, 'dev.off()'
    
    
	print >>ROUT, '\n'
	
	print >>ROUT, '%s(\"%s.%s\")' % (format.lower(),file_prefix + ".curves",format.lower())	
	print >>ROUT, "x=1:%d" % (colNum) 
	print >>ROUT, 'icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(%d)' % (len(names))
	
	if len(names) == 1:
		print >>ROUT, "plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0])
	
	elif  len(names) >=2 and len(names) <=6:  
		print >>ROUT, "plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0])
		for i in range(1,len(names)):
			print >>ROUT, "lines(x,%s,type='l',col=icolor[%d])" % (names[i], i+1)
		print >>ROUT, "legend(0,1,fill=icolor[%d:%d], legend=c(%s))" % (1,len(names), ','.join([ "'" + str(n) + "'" for n in names]))
	
	elif len(names) > 6:
		print >>ROUT, 'layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), 4, 4, byrow = TRUE))'
		print >>ROUT, "plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0])
		for i in range(1,len(names)):
			print >>ROUT, "lines(x,%s,type='l',col=icolor[%d])" % (names[i], i+1)
		print >>ROUT, 'par(mar=c(1,0,2,1))'
		print >>ROUT, 'plot.new()'
		print >>ROUT, "legend(0,1,fill=icolor[%d:%d],legend=c(%s))" % (1,len(names), ','.join([ "'" + str(n) + "'" for n in names]))

	print >>ROUT, 'dev.off()'
	ROUT.close()

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input file(s) in BAM format. "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files. 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam file (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools.')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. [required]")
	parser.add_option("-l","--minimum_length",action="store",type="int",default=100, dest="min_mRNA_length",help="Minimum mRNA length (bp). mRNA smaller than \"min_mRNA_length\" will be skipped. default=%default")
	parser.add_option("-f","--format",action="store",type="string",dest="output_format", default='pdf', help="Output file format, 'pdf', 'png' or 'jpeg'. default=%default")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print >>sys.stderr, '\n\n' + options.ref_gene_model + " does NOT exists" + '\n'
		#parser.print_help()
		sys.exit(0)
	if options.min_mRNA_length < 100:
		print >>sys.stderr, 'The number specified to "-l" cannot be smaller than 100.' + '\n'
		sys.exit(0)
		
	OUT1 = open(options.output_prefix  + ".geneBodyCoverage.txt"	,'w')
	print >>OUT1, "Percentile\t" + '\t'.join([str(i) for i in range(1,101)])
		
	printlog("Read BED file (reference gene model) ...")
	gene_percentiles = genebody_percentile(refbed = options.ref_gene_model, mRNA_len_cut = options.min_mRNA_length)
		
	printlog("Get BAM file(s) ...")
	bamfiles = getBamFiles.get_bam_files(options.input_files)
	for f in bamfiles:
		print >>sys.stderr, "\t" + f
	
	file_container = []
	for bamfile in bamfiles:
		printlog("Processing " + basename(bamfile) + ' ...')
		cvg = genebody_coverage(bamfile, gene_percentiles)
		if len(cvg) == 0:
			print >>sys.stderr, "\nCannot get coverage signal from " + basename(bamfile) + ' ! Skip'
			continue
		tmp = valid_name(basename(bamfile).replace('.bam',''))	# scrutinize R identifer
		if file_container.count(tmp) == 0:
			print >>OUT1, tmp + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)])
		else:
			print >>OUT1, tmp + '.' + str(file_container.count(tmp)) + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)])
		file_container.append(tmp)
	OUT1.close()
	
	
	dataset=[]
	for line in open(options.output_prefix  + ".geneBodyCoverage.txt",'r'):
		line = line.strip()
		if line.startswith("Percentile"):
			continue
		f = line.split()
		name = f[0]
		dat = [float(i) for i in  f[1:]]
		skewness = pearson_moment_coefficient(dat)
		dataset.append((name, [(i -min(dat))/(max(dat) - min(dat)) for i in dat], skewness))	
	dataset.sort(key = operator.itemgetter(2), reverse=True)
	
	print >>sys.stderr, "\n\n"
	print >>sys.stderr, "\tSample\tSkewness"
	for a,b,c in dataset:
		print >>sys.stderr, '\t' + a + '\t' + str(c)
	Rcode_write(dataset, options.output_prefix + '.geneBodyCoverage', format = options.output_format)
	
	printlog("Running R script ...")
	try:
		subprocess.call("Rscript " + options.output_prefix + '.geneBodyCoverage.r',shell=True)
	except:
		print >>sys.stderr, "Cannot generate pdf file from " + options.output_prefix + '.geneBodyCoverage.r'
		pass
	
	
if __name__ == '__main__':
	main()
	