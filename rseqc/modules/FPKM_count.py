#!/usr/bin/env python

'''
calculate raw read count, FPM (fragment per million) and FPKM (fragment per million mapped reads per kilobase exon) for each gene in BED file.
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

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
from qcmodule import bam_cigar
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

def build_range(refgene):
	'''build ranges for exonic region'''
	ranges={}
	for line in open(refgene,'r'):
		try:
			if line.startswith(('#','track','browser')):continue

			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0].upper()
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5].replace(" ","_")
				
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue	

		for st,end in zip(exon_starts,exon_ends):	
			if chrom not in ranges:
				ranges[chrom] = Intersecter()
			ranges[chrom].add_interval( Interval( st, end ) )
	return ranges
	

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM format (SAM is not supported). [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat. [required]")
	parser.add_option("-d","--strand",action="store",type="string",dest="strand_rule",default=None,help="How read(s) were stranded during sequencing. For example: --strand='1++,1--,2+-,2-+' means that this is a pair-end, strand-specific RNA-seq, and the strand rule is: read1 mapped to '+' => parental gene on '+'; read1 mapped to '-' => parental gene on '-'; read2 mapped to '+' => parental gene on '-'; read2 mapped to '-' => parental gene on '+'.  If you are not sure about the strand rule, run \'infer_experiment.py' default=%default (Not a strand specific RNA-seq data)")
	parser.add_option("-u","--skip-multi-hits",action="store_true",dest="skip_multi",help="How to deal with multiple hit reads. Presence this option renders program to skip multiple hits reads.")
	parser.add_option("-e","--only-exonic",action="store_true",dest="only_exon",help="How to count total reads. Presence of this option renders program only used exonic (UTR exons and CDS exons) reads, otherwise use all reads.")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")
	parser.add_option("-s","--single-read",action="store",type="float",dest="single_read",default=1,help="How to count read-pairs that only have one end mapped. 0: ignore it. 0.5: treat it as half fragment. 1: treat it as whole fragment. default=%default")
	
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file and options.refgene_bed):
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

	obj = SAM.ParseBAM(options.input_file)
	OUT = open(options.output_prefix + '.FPKM.xls','w')

	#++++++++++++++++++++++++++++++++++++determine strand rule
	strandRule={}
	if options.strand_rule is None:													# Not strand-specific
		pass																
	elif len(options.strand_rule.split(',')) ==4:									#PairEnd, strand-specific
		for i in options.strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
	elif len(options.strand_rule.split(',')) ==2:									#singeEnd, strand-specific
		for i in options.strand_rule.split(','):strandRule[i[0]]=i[1]
	else:
		print >>sys.stderr, "Unknown value of option :'strand_rule' " +  options.strand_rule
		sys.exit(1)	

	
	#++++++++++++++++++++++++++++++++++++counting fragments
	print >>sys.stderr, "Extract exon regions from  "+ options.refgene_bed + '...'
	gene_ranges = build_range( options.refgene_bed)
	print >>sys.stderr, "Counting total fragment ... ",
	
	total_frags =0.0
	exonic_frags = 0.0
	
	try:
		while(1):
			aligned_read = obj.samfile.next()
			if aligned_read.is_qcfail:continue			#skip low quanlity					
			if aligned_read.is_duplicate:continue		#skip duplicate read
			if aligned_read.is_secondary:continue		#skip non primary hit
			if options.skip_multi:
				if aligned_read.mapq < options.map_qual:
					continue
			try:
				chrom = obj.samfile.getrname(aligned_read.tid).upper()
			except:
				continue
			read_st = aligned_read.pos
			read_end = read_st + aligned_read.rlen 	#not exactly the end position in case of splicing, insertion,etc

			if not aligned_read.is_paired:	# if read is NOT paired in sequencing (single-end sequencing)
				total_frags += 1
				if (chrom in gene_ranges) and len(gene_ranges[chrom].find(read_st, read_end)) >0:
					exonic_frags += 1
			elif aligned_read.is_paired:	# for pair-end sequencing
				if aligned_read.is_read2: continue	# only count read1
				mate_st = aligned_read.pnext
				mate_end = mate_st + aligned_read.rlen
				
				if aligned_read.is_unmapped:	#read1 unmapped
					if aligned_read.mate_is_unmapped: continue	#both unmap
					else:	#read2 is mapped
						total_frags += options.single_read
						if (chrom in gene_ranges)  and (len(gene_ranges[chrom].find(mate_st, mate_end)) >0 ):
							exonic_frags += options.single_read
				else:
					if aligned_read.mate_is_unmapped:
						total_frags += options.single_read
						if (chrom in gene_ranges) and (len(gene_ranges[chrom].find(read_st, read_end)) >0 ):
							exonic_frags += options.single_read
					else:
						total_frags += 1
						if (chrom in gene_ranges) and (len(gene_ranges[chrom].find(read_st, read_end)) >0 ) and (len(gene_ranges[chrom].find(mate_st, mate_end)) >0 ):
							exonic_frags += 1
					
	except StopIteration:
		print >>sys.stderr, "Done"
	print >>sys.stderr, "Total fragment = %-20s" % (str(total_frags))
	print >>sys.stderr, "Total exonic fragment = %-20s" % (str(exonic_frags))

	if total_frags >0 and exonic_frags > 0:
		if options.only_exon:
			denominator = exonic_frags
		else:
			denominator = total_frags
	else:
		print >>sys.stderr, "Total tags cannot be 0 or negative number"
		sys.exit(1)
	
	
	#++++++++++++++++++++++++++++++++++++++++++++++++
	obj = SAM.ParseBAM(options.input_file)
	print >>OUT, '\t'.join(('#chrom','st','end','accession','mRNA_size','gene_strand','Frag_count','FPM','FPKM'))

	gene_finished=0

	#calculate raw count, FPM, FPKM for each gene
	for line in open(options.refgene_bed,'r'):
		frag_count_f = 0.0
		frag_count_r = 0.0
		frag_count_fr = 0.0
		mRNA_size = 0.0
		exon_ranges= Intersecter()
		if line.startswith(('#','track','browser')):continue   
		fields = line.split()
		chrom     = fields[0]
		tx_start  = int( fields[1] )
		tx_end    = int( fields[2] )
		geneName      = fields[3]
		gstrand    = fields[5].replace(" ","_")
	    	
		exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
		exon_starts = map((lambda x: x + tx_start ), exon_starts)
		exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
		exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
		for st,end in zip(exon_starts,exon_ends):
			mRNA_size += (end - st)
			exon_ranges.add_interval( Interval( st, end ) )
		
		
		# extract reads mapped gene region
		try:
			alignedReads = obj.samfile.fetch(chrom,tx_start,tx_end)
		except:
			continue
		for aligned_read in alignedReads:
			flag=0
			if aligned_read.is_qcfail:continue			#skip low quanlity					
			if aligned_read.is_duplicate:continue		#skip duplicate read
			if aligned_read.is_secondary:continue		#skip non primary hit		
			if options.skip_multi:
				if aligned_read.mapq < options.map_qual:
					continue
			
			#single end sequencing
			if not aligned_read.is_paired:
				frag_st = aligned_read.pos
				frag_end = read_st + aligned_read.rlen 	#not exactly the end position in case of splicing, insertion,etc
				if aligned_read.is_reverse:
					strand_key = '-'
				else:
					strand_key = '+'
				
				if len(exon_ranges.find(frag_st, frag_end)) > 0:
					if options.strand_rule is None:
						frag_count_fr += 1
					elif strand_key in strandRule and strandRule[strand_key] == '+':
						frag_count_f += 1
					elif strand_key in strandRule and strandRule[strand_key] == '-':
						frag_count_r += 1
			
			# pair-end sequencing
			if aligned_read.is_paired:
				frag_st = aligned_read.pos
				frag_end = aligned_read.pnext
				if len(exon_ranges.find(frag_st, frag_st +1 ))  < 1 and len(exon_ranges.find(frag_end, frag_end +1 )) < 1:
					continue
				if aligned_read.is_read2:
					continue
				if aligned_read.is_reverse:
					strand_key = '1-'
				else:
					strand_key = '1+'								
				
				if options.strand_rule is None:
					if aligned_read.is_unmapped:
						if aligned_read.mate_is_unmapped:	# both unmapped
							continue
						else:	#only read2 mapped
							frag_count_fr += options.single_read
					else:
						if aligned_read.mate_is_unmapped:	# only read1 mapped
							frag_count_fr += options.single_read
						else:	#both mapped
							frag_count_fr += 1
				else:
					if strand_key in strandRule and strandRule[strand_key] == '+':
						if aligned_read.is_unmapped:
							if aligned_read.mate_is_unmapped:	# both unmapped
								continue
							else:	#only read2 mapped
								frag_count_f += options.single_read
						else:
							if aligned_read.mate_is_unmapped:	# only read1 mapped
								frag_count_f += options.single_read
							else:	#both mapped
								frag_count_f += 1
					if strand_key in strandRule and strandRule[strand_key] == '-':
						if aligned_read.is_unmapped:
							if aligned_read.mate_is_unmapped:	# both unmapped
								continue
							else:	#only read2 mapped
								frag_count_r += options.single_read
						else:
							if aligned_read.mate_is_unmapped:	# only read1 mapped
								frag_count_r += options.single_read
							else:	#both mapped
								frag_count_r += 1
		
		FPM_fr =  frag_count_fr * 1000000 / denominator
		FPM_f =  frag_count_f * 1000000 / denominator
		FPM_r =  frag_count_r * 1000000 / denominator
		FPKM_fr =  frag_count_fr * 1000000000 / (denominator * mRNA_size)
		FPKM_f =  frag_count_f * 1000000000 / (denominator * mRNA_size)
		FPKM_r =  frag_count_r * 1000000000 / (denominator * mRNA_size)
		
		
		if options.strand_rule is None:
			print >>OUT, '\t'.join([str(i) for i in (chrom, tx_start, tx_end, geneName, mRNA_size, gstrand, frag_count_fr, FPM_fr,FPKM_fr)])
		else:
			if gstrand == '+':
				print >>OUT, '\t'.join([str(i) for i in (chrom, tx_start, tx_end, geneName, mRNA_size, gstrand, frag_count_f, FPM_f, FPKM_f )])
			elif gstrand == '-':
				print >>OUT, '\t'.join([str(i) for i in (chrom, tx_start, tx_end, geneName, mRNA_size, gstrand, frag_count_r, FPM_r, FPKM_r)])

		gene_finished +=1
		print >>sys.stderr, " %d transcripts finished\r" % (gene_finished),

if __name__ == '__main__':
	main()
