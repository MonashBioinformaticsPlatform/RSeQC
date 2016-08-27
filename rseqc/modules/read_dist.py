#!/usr/bin/env python
'''
Check reads distribution over exon, intron, UTR, intergenic ... etc
The following reads will be skipped:
	qc_failed
	PCR duplicate
	Unmapped
	Non-primary (or secondary)	
'''

#import built-in modules
import sys
import logging

# it looks like these imports aren't used in this module
#import logging
#import os
#import string
#import warnings
#import collections
#import math
#import sets

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *

#import my own modules
from rseqc.parsers import BED
from rseqc.qcmodule import SAM
from rseqc.qcmodule import bam_cigar
from rseqc.parsers.BedWrapper import BedWrapper

def cal_size(list):
	'''calcualte bed list total size'''
	size=0
	for l in list:
		size += l[2] - l[1]
	return size

def foundone(chrom,ranges, start, end):
	found = 0
	if chrom in ranges:
		found = len(ranges[chrom].find(start,end))
	return found

def build_bitsets(list):
	'''build intevalTree from list'''
	ranges={}
	for l in list:
		chrom =l[0].upper()
		start = int(l[1])
		end = int(l[2])
		if chrom not in ranges:
			ranges[chrom] = Intersecter()
		ranges[chrom].add_interval( Interval( start, end ) )
	return ranges

def process_gene_model(gene_model, file_type):
    logger = logging.getLogger(__name__)
    logger.info(' started procecssing annotation file > %s ...' % gene_model)

    gene_obj = None

    if file_type == 'gtf':
        gene_obj = BedWrapper(gene_model)
    if file_type == 'bed':
        gene_obj = BED.ParseBED(gene_model)

    # get certain aspects from that object
    # getUTR method retuns list of lists where nested list has 6 elements per gene
    utr_3 = gene_obj.getUTR(utr = 3)
    utr_5 = gene_obj.getUTR(utr = 5)
    cds_exon = gene_obj.getCDSExon()
    # getIntron method returns list of lists where nested list has 3 elements
    intron = gene_obj.getIntron()
    # I think the unionBed3 tries to collaps all overlapping regions... 
    intron = BED.unionBed3(intron)
    cds_exon = BED.unionBed3(cds_exon)
    utr_5 = BED.unionBed3(utr_5)
    utr_3 = BED.unionBed3(utr_3)
    
    utr_5 = BED.subtractBed3(utr_5, cds_exon)
    utr_3 = BED.subtractBed3(utr_3, cds_exon)
    intron = BED.subtractBed3(intron, cds_exon)
    intron = BED.subtractBed3(intron, utr_5)
    intron = BED.subtractBed3(intron, utr_3)
    
    intergenic_up_1kb = gene_obj.getIntergenic(direction="up",size=1000)
    intergenic_down_1kb = gene_obj.getIntergenic(direction="down",size=1000)
    intergenic_up_5kb = gene_obj.getIntergenic(direction="up",size=5000)
    intergenic_down_5kb = gene_obj.getIntergenic(direction="down",size=5000)	
    intergenic_up_10kb = gene_obj.getIntergenic(direction="up",size=10000)
    intergenic_down_10kb = gene_obj.getIntergenic(direction="down",size=10000)
    
    #merge integenic region
    intergenic_up_1kb=BED.unionBed3(intergenic_up_1kb)
    intergenic_up_5kb=BED.unionBed3(intergenic_up_5kb)
    intergenic_up_10kb=BED.unionBed3(intergenic_up_10kb)
    intergenic_down_1kb=BED.unionBed3(intergenic_down_1kb)
    intergenic_down_5kb=BED.unionBed3(intergenic_down_5kb)
    intergenic_down_10kb=BED.unionBed3(intergenic_down_10kb)	
    
    #purify intergenic region
    intergenic_up_1kb=BED.subtractBed3(intergenic_up_1kb,cds_exon)
    intergenic_up_1kb=BED.subtractBed3(intergenic_up_1kb,utr_5)
    intergenic_up_1kb=BED.subtractBed3(intergenic_up_1kb,utr_3)
    intergenic_up_1kb=BED.subtractBed3(intergenic_up_1kb,intron)
    intergenic_down_1kb=BED.subtractBed3(intergenic_down_1kb,cds_exon)
    intergenic_down_1kb=BED.subtractBed3(intergenic_down_1kb,utr_5)
    intergenic_down_1kb=BED.subtractBed3(intergenic_down_1kb,utr_3)
    intergenic_down_1kb=BED.subtractBed3(intergenic_down_1kb,intron)	
    
    #purify intergenic region
    intergenic_up_5kb=BED.subtractBed3(intergenic_up_5kb,cds_exon)
    intergenic_up_5kb=BED.subtractBed3(intergenic_up_5kb,utr_5)
    intergenic_up_5kb=BED.subtractBed3(intergenic_up_5kb,utr_3)
    intergenic_up_5kb=BED.subtractBed3(intergenic_up_5kb,intron)
    intergenic_down_5kb=BED.subtractBed3(intergenic_down_5kb,cds_exon)
    intergenic_down_5kb=BED.subtractBed3(intergenic_down_5kb,utr_5)
    intergenic_down_5kb=BED.subtractBed3(intergenic_down_5kb,utr_3)
    intergenic_down_5kb=BED.subtractBed3(intergenic_down_5kb,intron)	
    
    #purify intergenic region
    intergenic_up_10kb=BED.subtractBed3(intergenic_up_10kb,cds_exon)
    intergenic_up_10kb=BED.subtractBed3(intergenic_up_10kb,utr_5)
    intergenic_up_10kb=BED.subtractBed3(intergenic_up_10kb,utr_3)
    intergenic_up_10kb=BED.subtractBed3(intergenic_up_10kb,intron)
    intergenic_down_10kb=BED.subtractBed3(intergenic_down_10kb,cds_exon)
    intergenic_down_10kb=BED.subtractBed3(intergenic_down_10kb,utr_5)
    intergenic_down_10kb=BED.subtractBed3(intergenic_down_10kb,utr_3)
    intergenic_down_10kb=BED.subtractBed3(intergenic_down_10kb,intron)	
    
    #build intervalTree 
    cds_exon_ranges = build_bitsets(cds_exon)
    utr_5_ranges = build_bitsets(utr_5)
    utr_3_ranges = build_bitsets(utr_3)
    intron_ranges = build_bitsets(intron)
    
    interg_ranges_up_1kb_ranges = build_bitsets(intergenic_up_1kb)
    interg_ranges_up_5kb_ranges = build_bitsets(intergenic_up_5kb)
    interg_ranges_up_10kb_ranges = build_bitsets(intergenic_up_10kb)
    interg_ranges_down_1kb_ranges = build_bitsets(intergenic_down_1kb)
    interg_ranges_down_5kb_ranges = build_bitsets(intergenic_down_5kb)
    interg_ranges_down_10kb_ranges = build_bitsets(intergenic_down_10kb)
    
    exon_size = cal_size(cds_exon)
    intron_size = cal_size(intron)
    utr3_size = cal_size(utr_3)
    utr5_size = cal_size(utr_5)
    int_up1k_size = cal_size(intergenic_up_1kb)
    int_up5k_size = cal_size(intergenic_up_5kb)
    int_up10k_size = cal_size(intergenic_up_10kb)
    int_down1k_size = cal_size(intergenic_down_1kb)
    int_down5k_size = cal_size(intergenic_down_5kb)
    int_down10k_size = cal_size(intergenic_down_10kb)
    
    logger.info(' done processing annotation file')
    
    return (cds_exon_ranges,
            utr_5_ranges,
            utr_3_ranges,
            intron_ranges,
    	interg_ranges_up_1kb_ranges,
            interg_ranges_up_5kb_ranges,
            interg_ranges_up_10kb_ranges,
    	interg_ranges_down_1kb_ranges,
            interg_ranges_down_5kb_ranges,
            interg_ranges_down_10kb_ranges,
    	exon_size,
            intron_size,
            utr5_size,
            utr3_size,
    	int_up1k_size,
            int_up5k_size,
            int_up10k_size,
    	int_down1k_size,
            int_down5k_size,
            int_down10k_size
            )

def main(input_file, gene_models, file_type):

        logger = logging.getLogger(__name__)
        logger.info(" these are your input parameters: alignment file > %s | annotation file > %s | file type > %s" % (input_file, gene_models, file_type))

	# build bitset
	(cds_exon_r,
         intron_r,
         utr_5_r,
         utr_3_r,
	 intergenic_up_1kb_r,
         intergenic_up_5kb_r,
         intergenic_up_10kb_r,
	 intergenic_down_1kb_r,
         intergenic_down_5kb_r,
         intergenic_down_10kb_r,
	 cds_exon_base,
         intron_base,
         utr_5_base,
         utr_3_base,
	 intergenic_up1kb_base,
         intergenic_up5kb_base,
         intergenic_up10kb_base,
	 intergenic_down1kb_base,
         intergenic_down5kb_base,
         intergenic_down10kb_base) = process_gene_model(gene_models, file_type)
	
	intron_read=0
	cds_exon_read=0
	utr_5_read=0
	utr_3_read=0
	
	intergenic_up1kb_read=0
	intergenic_down1kb_read=0
	intergenic_up5kb_read=0
	intergenic_down5kb_read=0
	intergenic_up10kb_read=0
	intergenic_down10kb_read=0
		
	totalReads=0
	totalFrags=0
	unAssignFrags=0

	obj = SAM.ParseBAM(input_file)
	
	R_qc_fail=0
	R_duplicate=0
	R_nonprimary=0
	R_unmap=0
	
        logger.info(' started processing alignemnt file > %s ..' % input_file)

	try:
		while(1):
			aligned_read = obj.samfile.next()
			if aligned_read.is_qcfail:			#skip QC fail read
				R_qc_fail +=1
				continue
			if aligned_read.is_duplicate:		#skip duplicate read
				R_duplicate +=1
				continue
			if aligned_read.is_secondary:		#skip non primary hit
				R_nonprimary +=1
				continue
			if aligned_read.is_unmapped:		#skip unmap read
				R_unmap +=1
				continue		
			totalReads +=1
			chrom = obj.samfile.getrname(aligned_read.tid)
			chrom = chrom.upper()
			exons = bam_cigar.fetch_exon(chrom, aligned_read.pos, aligned_read.cigar)
			totalFrags += len(exons)
			
			for exn in exons:
				#print chrom + '\t' + str(exn[1]) + '\t' + str(exn[2])
				mid = int(exn[1]) + int((int(exn[2]) - int(exn[1]))/2)
				if foundone(chrom,cds_exon_r,mid,mid) > 0:
					cds_exon_read += 1
					continue
				elif foundone(chrom,utr_5_r,mid,mid) >0 and foundone(chrom,utr_3_r,mid,mid) == 0:
					utr_5_read += 1
					continue
				elif foundone(chrom,utr_3_r,mid,mid) >0 and foundone(chrom,utr_5_r,mid,mid) == 0:
					utr_3_read += 1
					continue
				elif foundone(chrom,utr_3_r,mid,mid) >0 and foundone(chrom,utr_5_r,mid,mid) > 0:
					unAssignFrags +=1
					continue
				elif foundone(chrom,intron_r,mid,mid) > 0:
					intron_read += 1
					continue
				elif foundone(chrom,intergenic_up_10kb_r,mid,mid) >0 and foundone(chrom,intergenic_down_10kb_r,mid,mid) > 0:
					unAssignFrags +=1
					continue					
				elif foundone(chrom,intergenic_up_1kb_r,mid,mid) >0:
					intergenic_up1kb_read += 1
					intergenic_up5kb_read += 1
					intergenic_up10kb_read += 1
				elif foundone(chrom,intergenic_up_5kb_r,mid,mid) >0:
					intergenic_up5kb_read += 1
					intergenic_up10kb_read += 1
				elif foundone(chrom,intergenic_up_10kb_r,mid,mid) >0:
					intergenic_up10kb_read += 1
				
				elif foundone(chrom,intergenic_down_1kb_r,mid,mid) >0:
					intergenic_down1kb_read += 1
					intergenic_down5kb_read += 1
					intergenic_down10kb_read += 1
				elif foundone(chrom,intergenic_down_5kb_r,mid,mid) >0:
					intergenic_down5kb_read += 1
					intergenic_down10kb_read += 1
				elif foundone(chrom,intergenic_down_10kb_r,mid,mid) >0:
					intergenic_down10kb_read += 1	
				else:
					unAssignFrags +=1
	except StopIteration:
		logger.info(' done processing alignment file')

	print "%-30s%d" % ("Total Reads",totalReads)
	print  "%-30s%d" % ("Total Tags",totalFrags)
	print  "%-30s%d" % ("Total Assigned Tags",totalFrags-unAssignFrags)
	
	print  "====================================================================="
	print  "%-20s%-20s%-20s%-20s" % ('Group','Total_bases','Tag_count','Tags/Kb')
	print  "%-20s%-20d%-20d%-18.2f" % ('CDS_Exons',cds_exon_base,cds_exon_read,cds_exon_read*1000.0/(cds_exon_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("5'UTR_Exons",utr_5_base,utr_5_read, utr_5_read*1000.0/(utr_5_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("3'UTR_Exons",utr_3_base,utr_3_read, utr_3_read*1000.0/(utr_3_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("Introns",intron_base,intron_read,intron_read*1000.0/(intron_base+1))
	
	print  "%-20s%-20d%-20d%-18.2f" % ("TSS_up_1kb",intergenic_up1kb_base, intergenic_up1kb_read, intergenic_up1kb_read*1000.0/(intergenic_up1kb_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("TSS_up_5kb",intergenic_up5kb_base, intergenic_up5kb_read, intergenic_up5kb_read*1000.0/(intergenic_up5kb_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("TSS_up_10kb",intergenic_up10kb_base, intergenic_up10kb_read, intergenic_up10kb_read*1000.0/(intergenic_up10kb_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("TES_down_1kb",intergenic_down1kb_base, intergenic_down1kb_read, intergenic_down1kb_read*1000.0/(intergenic_down1kb_base+1))
	print  "%-20s%-20d%-20d%-18.2f" % ("TES_down_5kb",intergenic_down5kb_base, intergenic_down5kb_read, intergenic_down5kb_read*1000.0/(intergenic_down5kb_base+1))	
	print  "%-20s%-20d%-20d%-18.2f" % ("TES_down_10kb",intergenic_down10kb_base, intergenic_down10kb_read, intergenic_down10kb_read*1000.0/(intergenic_down10kb_base+1))
	print  "====================================================================="
