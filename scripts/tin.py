#!/usr/bin/env python
'''-------------------------------------------------------------------------------------------------
This program calculates transcript integrity number (TIN) for each transcript (or gene) in
BED file. TIN is conceptually similar to RIN (RNA integrity number) but provides transcript
level measurement of RNA quality and is more sensitive to measure low quality RNA samples:

1) TIN score of a transcript is used to measure the RNA integrity of the transcript.
2) Median TIN score across all transcripts can be used to measure RNA integrity  of that
   "RNA sample".
3) TIN ranges from 0 (the worst) to 100 (the best). TIN = 60 means: 60% of the transcript
   has been covered if the reads coverage were uniform.
4) TIN will be assigned to 0 if the transcript has no coverage or covered reads is fewer than
   cutoff.
-------------------------------------------------------------------------------------------------'''
import sys,os
import math,random
from optparse import OptionParser
from qcmodule import SAM
from qcmodule import getBamFiles
from numpy import mean,median,std
from time import strftime
import pysam
from bx.intervals import *

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
	'''
	print mesg into stderr with time string appending to it.
	'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	print >>sys.stderr,mesg

def uniqify(seq):
	'''
	duplicated members only keep one copy. [1,2,2,3,3,4] => [1,2,3,4].
	'''
	seen = set()
	return [x for x in seq if x not in seen and not seen.add(x)]

def shannon_entropy(arg):
	'''
	calculate shannon's H = -sum(P*log(P)). arg is a list of float numbers. Note we used
	natural log here.
	'''
	lst_sum = sum(arg)
	entropy = 0.0
	for i in arg:
		entropy += (i/lst_sum) * math.log(i/lst_sum)
	if entropy == 0:
		return 0
	else:
		return -entropy

def build_bitsets(list):
	'''
	build intevalTree from list
	'''
	ranges={}
	for l in list:
		chrom =l[0]
		st = l[1]
		end = l[2]
		if chrom not in ranges:
			ranges[chrom] = Intersecter()
		ranges[chrom].add_interval( Interval( st, end ) )
	return ranges

def union_exons(refbed):
	'''
	take the union of all exons defined in refbed file and build bitset
	'''
	from qcmodule import BED
	tmp = BED.ParseBED(refbed)
	all_exons = tmp.getExon()
	unioned_exons = BED.unionBed3(all_exons)
	exon_ranges = build_bitsets(unioned_exons)
	return exon_ranges

def estimate_bg_noise(chrom, tx_st, tx_end, samfile, e_ranges):
	'''
	estimate background noise level for a particular transcript
	'''
	intron_sig = 0.0	# reads_num * reads_len
	alignedReads = samfile.fetch(chrom,tx_st,tx_end)
	for aligned_read in alignedReads:
		if aligned_read.is_qcfail:continue 
		if aligned_read.is_unmapped:continue
		if aligned_read.is_secondary:continue
		read_start = aligned_read.pos
		if read_start < tx_st: continue
		if read_start >= tx_end: continue
		read_len = aligned_read.qlen
		if len(e_ranges[chrom].find(read_start, read_start + read_len)) > 0:
			continue
		intron_sig += read_len
	return intron_sig

def genomic_positions(refbed, sample_size):
	''' 
	return genomic positions of each nucleotide in mRNA. sample_size: Number of nucleotide
	positions sampled from mRNA.
	'''
	if refbed is None:
		print >>sys.stderr,"You must specify a bed file representing gene model\n"
		exit(0)
	
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
			cdsStart = int(fields[6]) + 1	#convert to 1-based
			cdsEnd = int(fields[7])
			exon_count = int(fields[9])
			mRNA_size = sum([int(i) for i in fields[10].strip(',').split(',')])
			geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
				
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
			intron_size = tx_end - tx_start - mRNA_size
			if intron_size < 0:
				intron_size = 0
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue
		
		chose_bases=[tx_start+1, tx_end]
		exon_bounds = []
		gene_all_base = []
		if mRNA_size <= sample_size:	# return all bases of mRNA
			for st,end in zip(exon_starts,exon_ends):
				chose_bases.extend(range(st+1,end+1))		#1-based coordinates on genome, include exon boundaries
			yield (geneName, chrom, tx_start, tx_end, intron_size, chose_bases)
		elif mRNA_size > sample_size:
			step_size = int(mRNA_size/sample_size)
			for st,end in zip(exon_starts,exon_ends):
				gene_all_base.extend(range(st+1,end+1))
				exon_bounds.append(st+1)
				exon_bounds.append(end)
			indx = range(0, len(gene_all_base), step_size)
			chose_bases = [gene_all_base[i] for i in indx]
			yield (geneName, chrom, tx_start, tx_end, intron_size, uniqify(exon_bounds + chose_bases))
		
def check_min_reads(samfile,chrom,tx_st, tx_end, cutoff):
	'''
	make sure the gene has minimum reads coverage. if cutoff = 10, each gene must have
	10 *different* reads.
	'''
	tmp = False
	read_count = set()
	try:
		alignedReads = samfile.fetch(chrom,tx_st,tx_end)
		for aligned_read in alignedReads:
			if aligned_read.is_qcfail:continue 
			if aligned_read.is_unmapped:continue
			if aligned_read.is_secondary:continue
			read_start = aligned_read.pos
			if read_start < tx_st: continue
			if read_start >= tx_end: continue
			read_count.add(read_start)
			if len(read_count) > cutoff:	#no need to loop anymore
				tmp = True
				break
		return tmp
	except:
		return False

	
def genebody_coverage(samfile, chrom, positions, bg_level = 0):
	'''
	calculate coverage for each nucleotide in *positions*. some times len(cvg) < len(positions)
	because positions where there is no mapped reads were ignored.
	'''
	cvg = []
	start = positions[0] - 1
	end = positions[-1]
	
	try:
		for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions: continue
			if pileupcolumn.n == 0:
				cvg.append(0.0)
				continue
			cover_read = 0.0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				#if pileupread.alignment.is_duplicate:continue
				cover_read +=1.0
			cvg.append(cover_read)
	except:
		cvg = []
	
	if bg_level <= 0:
		return cvg
	else:
		tmp = []
		for i in cvg:
			subtracted_sig = int(i - bg_level)
			if subtracted_sig > 0:
				tmp.append(subtracted_sig)
			else:
				tmp.append(0)
		return tmp


def tin_score(cvg, l):
	'''calcualte TIN score'''
	
	if len(cvg) == 0:
		tin = 0
		return tin
	
	cvg_eff = [float(i) for i in cvg if float(i) > 0]	#remove positions with 0 read coverage
	entropy = shannon_entropy(cvg_eff)
	
	tin = 100*(math.exp(entropy)) / l
	return tin


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input BAM file(s). "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files (no spaces allowed). 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam files (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools. [required]')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]")
	parser.add_option("-c","--minCov",action="store",type="int",dest="minimum_coverage",default=10,help="Minimum number of read mapped to a transcript. default=%default")
	parser.add_option("-n","--sample-size",action="store",type="int",dest="sample_size",default=100,help="Number of equal-spaced nucleotide positions picked from mRNA. Note: if this number is larger than the length of mRNA (L), it will be halved until it's smaller than L. default=%default")
	parser.add_option("-s","--subtract-background",action="store_true",dest="subtract_bg",help="Subtract background noise (estimated from intronic reads). Only use this option if there are substantial intronic reads.")
	(options,args)=parser.parse_args()
	
	# if '-s' was set
	if options.subtract_bg:
		exon_ranges = union_exons(options.ref_gene_model)
		
	if options.sample_size < 0:
		print >>sys.stderr, "Number of nucleotide can't be negative"
		sys.exit(0)
	elif options.sample_size >1000:
		print >>sys.stderr, "Warning: '-n' is too large! Please try smaller '-n' valeu if program is running slow."
		
	if not (options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print >>sys.stderr, '\n\n' + options.ref_gene_model + " does NOT exists" + '\n'
		parser.print_help()
		sys.exit(0)
		
	printlog("Get BAM file(s) ...")
	bamfiles = sorted(getBamFiles.get_bam_files(options.input_files))
	
	if len(bamfiles) <= 0:
		print >>sys.stderr, "No BAM file found, exit."
		sys.exit(0)
	else:
		print >>sys.stderr, "Total %d BAM file(s):" % len(bamfiles)
		for f in bamfiles:
			print >>sys.stderr, "\t" + f	
	
	
	for f in bamfiles:
		printlog("Processing " + f)
		
		SUM = open(os.path.basename(f).replace('bam','') + 'summary.txt','w')
		print >>SUM, "\t".join(['Bam_file','TIN(mean)', 'TIN(median)','TIN(stdev)'])
		
		OUT = open(os.path.basename(f).replace('bam','') + 'tin.xls','w')
		print >>OUT, "\t".join(["geneID","chrom", "tx_start", "tx_end","TIN"])
		
		samfile = pysam.Samfile(f, "rb")
		sample_TINs = []	#sample level TIN, values are from different genes
		finish = 0
		noise_level = 0.0
		for gname, i_chr, i_tx_start, i_tx_end, intron_size, pick_positions in genomic_positions(refbed = options.ref_gene_model, sample_size = options.sample_size):	
			finish += 1
			
			# check minimum reads coverage
			if check_min_reads(samfile,i_chr,i_tx_start,i_tx_end,options.minimum_coverage) is not True:
				print >>OUT, '\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, 0.0)])
				continue
				
			# estimate background noise if '-s' was specified
			if options.subtract_bg:
				intron_signals = estimate_bg_noise(i_chr, i_tx_start, i_tx_end, samfile, exon_ranges)
				if intron_size > 0:
					noise_level = intron_signals/intron_size					

			coverage = genebody_coverage(samfile, i_chr,sorted(pick_positions), noise_level)
			
			#for a,b in zip(sorted(pick_positions),coverage):
			#	print str(a) + '\t' + str(b)
			
			tin1 = tin_score(cvg = coverage, l = len(pick_positions))
			sample_TINs.append(tin1)
			print >>OUT, '\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, tin1)])
			print >>sys.stderr, " %d transcripts finished\r" % (finish),
		
		print >>SUM, "\t".join( [str(i) for i in (os.path.basename(f), mean(sample_TINs), median(sample_TINs), std(sample_TINs))])
		OUT.close()
		SUM.close()
		samfile.close()

if __name__ == '__main__':
	
	main()

	
