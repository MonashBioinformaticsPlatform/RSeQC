#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
from rseqc.parsers.GTF import GeneModels
#from bx.bitset import *
#from bx.bitset_builders import *
#from bx.intervals import *

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
