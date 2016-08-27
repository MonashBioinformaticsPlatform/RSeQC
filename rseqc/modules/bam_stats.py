#!/usr/bin/env python
'''
Summarizing mapping statistics of a BAM or SAM file. 
'''

#import my own modules
from rseqc.qcmodule import SAM

__author__ = "Liguo Wang"
__copyright__ = "Copyleft."
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def main(input_file, mapq):
    obj = SAM.ParseBAM(input_file)
    obj.stat(q_cut = mapq)
