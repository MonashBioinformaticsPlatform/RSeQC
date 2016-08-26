#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
from rseqc.parsers.GTF import GeneModels
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

class BedWrapper(GeneModels):

    def getCDSExon(self):
        cds = []

        for c in self.get_cds():
            cds.append([c['chr'], c['start'], c['end']])

        return cds

    def getIntron(self):

        introns = []

        for i in self.get_introns():
            introns.append([i['chr'], i['start'], i['end']])

        return introns

    def getUTR(self, utr = 35):

        tprime_utrs = []

        if utr == 3: 
            for i in self.get_utrs(biotype = "protein_coding", tsl = 1, utr = 3):
                tprime_utrs.append([i['chr'], i['start'], i['end'], i['name'], '0', i['strand']])

        if utr == 5: 
            for i in self.get_utrs(biotype = "protein_coding", tsl = 1, utr = 5):
                tprime_utrs.append([i['chr'], i['start'], i['end'], i['name'], '0', i['strand']])

        return tprime_utrs

    def getIntergenic(self, direction = 'up', size = 1000):
        '''get intergenic regions. direction = up or down or both.'''
        
        regions = []
    
        for t in self.get_cds():
    
            chrom = t['chr']
            tx_start = t['start']
            tx_end = t['end']
            strand = t['strand']
            
            if(direction == "up" or direction == "both"):
                if strand == '-':
                    region_st = tx_end
                    region_end = tx_end + size
                else:
                    region_st = max(tx_start - size, 0)
                    region_end = tx_start
            
                regions.append([chrom, region_st, region_end])
            
            if (direction == "down" or direction == "both"):
                if strand == '-':
                    region_st = max(0,tx_start-size)
                    region_end = tx_start
                else:
                    region_st = tx_end
                    region_end = tx_end+size
                regions.append([chrom, region_st, region_end])
            
        return regions


class Helpers(object):

    def unionBed3(self, lst):
        '''Take the union of 3 column bed files. return a new list'''

        assert isinstance(lst, list)

        bitsets = binned_bitsets_from_list(lst)

        ret_lst = []
        for chrom in bitsets:
            bits = bitsets[chrom]
            end = 0
            while 1:
                start = bits.next_set( end )
                if start == bits.size:
                    break
                end = bits.next_clear( start )
                ret_lst.append([chrom, start, end])
        bitsets = dict()
        return ret_lst
    
    def subtractBed3(self, lst1, lst2):
        '''subtrack lst2 from lst1'''
        bitsets1 = binned_bitsets_from_list(lst1)
        bitsets2 = binned_bitsets_from_list(lst2)
        
        ret_lst=[]
        for chrom in bitsets1:  
            if chrom not in bitsets1:
                continue
            bits1 = bitsets1[chrom]
            if chrom in bitsets2:
                bits2 = bitsets2[chrom]
                bits2.invert()
                bits1.iand( bits2 )
            end=0
            while 1:
                start = bits1.next_set( end )
                if start == bits1.size:
                    break
                end = bits1.next_clear( start )
    
                ret_lst.append([chrom, start, end])
    
        bitsets1 = dict()
        bitsets2 = dict()
    
        return ret_lst

