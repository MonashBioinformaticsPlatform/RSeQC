#!/usr/bin/env python
'''manipulate BAM/SAM file.'''
#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
import random

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
import pysam
from rseqc.qcmodule import mystat
from rseqc.qcmodule import fasta
from rseqc.qcmodule import bam_cigar
from rseqc.qcmodule import BED
#changes to the paths

#changing history to this module
#05/26/2011: suppport multiple spliced mapped reads
#10/13/2011: saturation test for RNAs-eq data

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2013, Liguo Wang"
__credits__ = []
__license__ = "GPL"
__version__ = "2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"


class ParseSAM(object):
        '''This class provides fuctions to parsing/processing/transforming SAM format file
        Format of SAM file see: http://samtools.sourceforge.net/SAM1.pdf'''
        _reExpr1=re.compile(r'\s+')                             #>=1 spaces
        _reExpr2=re.compile(r'^\s*$')                   #blank line
        _splicedHit_pat = re.compile(r'(\d+)[M|N]',re.IGNORECASE)       #regular expression for spliced mapped reads
        _monoHit_pat = re.compile(r'^(\d+)M$',re.IGNORECASE)                            #regular expresion for Non-spliced mapped reads
        _insertionHit_pat =re.compile(r'\d+I',re.IGNORECASE)
        _deletionHit_pat =re.compile(r'\d+D',re.IGNORECASE)
        _softClipHi_pat =re.compile(r'\d+S',re.IGNORECASE)
        _hardClipHi_pat =re.compile(r'\d+H',re.IGNORECASE)
        _padHit_pat =re.compile(r'\d+P',re.IGNORECASE)
        _uniqueHit_pat = re.compile(r'[NI]H:i:1\b')

        def __init__(self,samFile):
                '''constructor'''
                if samFile == '-':
                        self.fileName = "STDIN"
                        self.f = sys.stdin
                else:
                        self.fileName = os.path.basename(samFile)
                        self.f=open(samFile,'r')

        def stat (self):
                '''Calculate mapping statistics'''              
                total_read=0
                
                pcr_duplicate =0
                low_qual =0
                secondary_hit =0
                
                unmapped_read1=0
                mapped_read1=0
                reverse_read1=0
                forward_read1=0
                
                unmapped_read2=0
                mapped_read2=0
                reverse_read2=0
                forward_read2=0
                
                _numSplitHit =0
                _numMonoHit =0
                _numInsertion =0
                _numDeletion =0

                minus_minus=0
                minus_plus =0
                plus_minus=0
                plus_plus=0
                paired=True
                
                unmap_SE=0
                map_SE=0
                reverse_SE=0
                forward_SE=0
                for line in self.f:
                        line=line.rstrip()
                        if line.startswith('@'):continue                                #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        
                        total_read +=1
                        field = line.split()
                        flagCode=string.atoi(field[1])
                        if  (flagCode & 0x0400 !=0):                                                    #PCR or optical duplicate
                                pcr_duplicate +=1
                                continue
                        if  (flagCode & 0x0200 !=0):                                                    #Low quality
                                low_qual +=1
                                continue
                        if  (flagCode & 0x0200 !=0):                                                    #Not primary alignment
                                secondary_hit +=1
                                continue
                        if (len(ParseSAM._splicedHit_pat.findall(field[5]))>1):_numSplitHit +=1                 #Splicing mapped reads                                                  
                        if (len(ParseSAM._splicedHit_pat.findall(field[5]))==1):_numMonoHit +=1                 #mono mapped reads                                      
                        if (ParseSAM._insertionHit_pat.search(field[5])):_numInsertion +=1                              #insertion in reads
                        if (ParseSAM._deletionHit_pat.search(field[5])):_numDeletion +=1                                #deletion in reads
                        


                        if (flagCode & 0x0001 !=0):                                                             #This is paired end sequencing
                                if (flagCode & 0x0040 != 0):                                            #1st read
                                        if (flagCode & 0x0004 != 0):unmapped_read1 +=1
                                        if (flagCode & 0x0004 == 0):mapped_read1 +=1
                                        if (flagCode & 0x0010 != 0):reverse_read1 +=1
                                        if (flagCode & 0x0010 == 0):forward_read1 +=1
                                                
                                if (flagCode & 0x0080 != 0):                                            #2nd read
                                        if (flagCode & 0x0004 != 0):unmapped_read2 +=1
                                        if (flagCode & 0x0004 == 0):mapped_read2 +=1
                                        if (flagCode & 0x0010 != 0):reverse_read2 +=1
                                        if (flagCode & 0x0010 == 0):forward_read2 +=1                                           
                                if      (flagCode & 0x0010 != 0 and flagCode & 0x0020 != 0):
                                        minus_minus +=1
                                if      (flagCode & 0x0010 != 0 and flagCode & 0x0020 == 0):
                                        minus_plus +=1
                                if      (flagCode & 0x0010 == 0 and flagCode & 0x0020 != 0):
                                        plus_minus +=1
                                if      (flagCode & 0x0010 == 0 and flagCode & 0x0020 == 0):
                                        plus_plus +=1                           
                        if (flagCode & 0x0001 ==0):                                                             #This is single end sequencing
                                paired=False
                                if (flagCode & 0x0004 != 0):
                                        unmap_SE +=1
                                if (flagCode & 0x0004 == 0):
                                        map_SE +=1
                                if (flagCode & 0x0010 != 0):
                                        reverse_SE +=1
                                if (flagCode & 0x0010 == 0):
                                        forward_SE +=1
                                
                if paired:                              
                        print("\n#==================================================", file=sys.stderr)
                        print("#================Report (pair-end)=================", file=sys.stderr)
                        print("%-25s%d" % ("Total Reads:",total_read), file=sys.stderr)
                        print("%-25s%d" % ("Total Mapped Reads:", (mapped_read1 + mapped_read2)), file=sys.stderr)
                        print("%-25s%d" % ("Total Unmapped Reads:",(unmapped_read1 + unmapped_read2)), file=sys.stderr)
                        print("%-25s%d" % ("PCR duplicate:",pcr_duplicate), file=sys.stderr)
                        print("%-25s%d" % ("QC-failed:",low_qual), file=sys.stderr)
                        print("%-25s%d" % ("Not primary mapping:",secondary_hit), file=sys.stderr)
                        print("\n", end=' ', file=sys.stderr)
                        print("%-25s%d" % ("Unmapped Read-1:",unmapped_read1), file=sys.stderr)
                        print("%-25s%d" % ("Mapped Read-1:",mapped_read1), file=sys.stderr)
                        print("%-25s%d" % ("  Forward (+):",forward_read1), file=sys.stderr)
                        print("%-25s%d" % ("  Reverse (-):",reverse_read1), file=sys.stderr)
                        
                        print("\n", end=' ', file=sys.stderr)
                        print("%-25s%d" % ("Unmapped Read-2:",unmapped_read2), file=sys.stderr)
                        print("%-25s%d" % ("Mapped Read-2:",mapped_read2), file=sys.stderr)
                        print("%-25s%d" % ("  Forward (+):",forward_read2), file=sys.stderr)
                        print("%-25s%d" % ("  Reverse (-):",reverse_read2), file=sys.stderr)
                        
                        print("\n", end=' ', file=sys.stderr)
                        print("%-25s%d" % ("Mapped to (+/-):",plus_minus), file=sys.stderr)
                        print("%-25s%d" % ("Mapped to (-/+):",minus_plus), file=sys.stderr)
                        print("%-25s%d" % ("Mapped to (+/+):",plus_plus), file=sys.stderr)
                        print("%-25s%d" % ("Mapped to (-/-):",minus_minus), file=sys.stderr)
                        print("\n", end=' ', file=sys.stderr)
                        print("%-25s%d" % ("Spliced Hits:",_numSplitHit), file=sys.stderr)
                        print("%-25s%d" % ("Non-spliced Hits:",_numMonoHit), file=sys.stderr)
                        print("%-25s%d" % ("Reads have insertion:",_numInsertion), file=sys.stderr)
                        print("%-25s%d" % ("Reads have deletion:",_numDeletion), file=sys.stderr)
                else:
                        print("\n#====================================================", file=sys.stderr)
                        print("#================Report (single-end)=================", file=sys.stderr)
                        print("%-25s%d" % ("Total Reads:",total_read), file=sys.stderr)
                        print("%-25s%d" % ("Total Mapped Reads:", map_SE), file=sys.stderr)
                        print("%-25s%d" % ("Total Unmapped Reads:",unmap_SE), file=sys.stderr)
                        print("%-25s%d" % ("PCR duplicate:",pcr_duplicate), file=sys.stderr)
                        print("%-25s%d" % ("QC-failed:",low_qual), file=sys.stderr)
                        print("%-25s%d" % ("Not primary mapping:",secondary_hit), file=sys.stderr)
                        print("%-25s%d" % ("froward (+):",forward_SE), file=sys.stderr)
                        print("%-25s%d" % ("reverse (-):",reverse_SE), file=sys.stderr)
                        print("\n", end=' ', file=sys.stderr)
                        print("%-25s%d" % ("Spliced Hits:",_numSplitHit), file=sys.stderr)
                        print("%-25s%d" % ("Non-spliced Hits:",_numMonoHit), file=sys.stderr)
                        print("%-25s%d" % ("Reads have insertion:",_numInsertion), file=sys.stderr)
                        print("%-25s%d" % ("Reads have deletion:",_numDeletion), file=sys.stderr)                       
                        
        def samTobed(self,outfile=None,mergePE=False):
                """Convert SAM file to BED file. BED file will be saved as xxx.sam.bed unless otherwise specified.
                 If mergePE=False, each read will be one bed entry. If mergePE=True, pair-end (if there are and on
                 the same chr) reads will be displayed in bed entry."""
                if outfile is None:
                        outfile=self.fileName + ".bed"

                print("\tWriting bed entries to\"",outfile,"\"...", end=' ', file=sys.stderr)
                FO=open(outfile,'w')
                for line in self.f:
                        if line.startswith(('@','track')):continue              #skip head lines
                        if ParseSAM._reExpr2.match(line):continue       #skip blank line
                        field=line.rstrip().split()
                        if (string.atoi(field[1]) & 0x0004)!=0: continue        #skip unmapped line
                        if (string.atoi(field[1]) & 0x0040)!=0:
                                mate="/1"
                        else:
                                mate="/2"
                        flag = string.atoi(field[1])
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]       #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        chrom = field[2]
                        chromStart = string.atoi(field[3])-1
                        chromEnd=chromStart
                        for i in comb:
                                chromEnd += i           
                        name = field[0] + mate
                        score = field[4]
                        if(flag & 0x0010)==0:
                                strand = '+'
                        else:
                                strand = '-'
                        thickStart = chromStart
                        thickEnd = chromEnd
                        itemRgb = "0,255,0"
                        blockCount = (len(comb) +1) /2
                        blockSize = []
                        for i in range(0,len(comb),2):
                                blockSize.append(str(comb[i]))
                        blockSizes = ','.join(blockSize)
                        blockStart=[]
                        for i in range(0,len(comb),2):
                                blockStart.append(str(sum(comb[:i])))
                        blockStarts = ','.join(blockStart)
                        print(string.join((str(i) for i in [chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts]),sep="\t"), file=FO)                     
                print("Done", file=sys.stderr)
                FO.close()
                self.f.seek(0)

                if mergePE:
                #creat another bed file. pair-end reads will be merged into single bed entry
                        print("Writing consoidated bed file ...", end=' ', file=sys.stderr)
                        bedfile = open(outfile,'r')
                        outfile_2 = outfile + ".consolidate.bed"
                        outfile_3 = outfile + '.filter'
                        FO = open(outfile_2,'w')
                        FOF = open(outfile_3,'w')
                        count={}
                        chr=collections.defaultdict(set)
                        txSt={}
                        txEnd={}
                        strand=collections.defaultdict(list)
                        blocks={}
                        sizes=collections.defaultdict(list)
                        starts=collections.defaultdict(list)
                        for line in bedfile:
                                field=line.strip().split()
                                field[3]=field[3].replace('/1','')
                                field[3]=field[3].replace('/2','')
                                if field[3] not in count:
                                        count[field[3]] = 1
                                        chr[field[3]].add(field[0])
                                        txSt[field[3]] = string.atoi(field[1])
                                        txEnd[field[3]] = string.atoi(field[2])
                                        strand[field[3]].append(field[5])
                                        blocks[field[3]] = string.atoi(field[9])
                                        sizes[field[3]].extend( field[10].split(',') )
                                        starts[field[3]].extend([string.atoi(i) + string.atoi(field[1]) for i in field[11].split(',') ])
                                else:
                                        count[field[3]] += 1
                                        chr[field[3]].add(field[0])
                                        if string.atoi(field[1]) < txSt[field[3]]:
                                                txSt[field[3]] = string.atoi(field[1]) 
                                        if string.atoi(field[2]) > txEnd[field[3]]:
                                                txEnd[field[3]] = string.atoi(field[2]) 
                                        blocks[field[3]] += string.atoi(field[9])
                                        strand[field[3]].append(field[5])
                                        sizes[field[3]].extend( field[10].split(',') )
                                        starts[field[3]].extend([string.atoi(i) + string.atoi(field[1]) for i in field[11].split(',') ])                                
                        #befile.seek(0)
                
                        for key in count:
                                st=[]   #real sorted starts
                                sz=[]   #real sorted sizes
                                if(count[key] ==1):                     #single-end read
                                        if(blocks[key] ==1):    #single end, single hit
                                                st = [i - txSt[key] for i in starts[key]]
                                                st = string.join([str(i) for i in st],',')
                                                print(chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key,"\t","11\t",strand[key][0],"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sizes[key],','),"\t",st, file=FO)
                                        else:
                                                st = [i - txSt[key] for i in starts[key]]       #single end, spliced hit
                                                st = string.join([str(i) for i in st],',')
                                                print(chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key,"\t","12\t",strand[key][0],"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sizes[key],','),"\t",st, file=FO)                                          
                                                                
                                elif(count[key]==2):    #pair-end read
                                        direction = string.join(strand[key],'/')        
                                        for i,j in sorted (zip(starts[key],sizes[key])):
                                                st.append(i-txSt[key])
                                                sz.append(j)
                                        #st=[string.atoi(i) for i in st]
                                        if(len(chr[key])==1):   #pair-end reads mapped to same chromosome
                                                if blocks[key] ==2:     #pair end, single hits
                                                        print(chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key + "|strand=" + direction + "|chrom=same","\t","21\t",'.',"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sz,','),"\t",string.join([str(i) for i in st],','), file=FO)
                                                elif blocks[key] >2:    #
                                                        print(chr[key].pop(),"\t",txSt[key],"\t",txEnd[key],"\t",key + "|strand=" + direction + "|chrom=same","\t","22\t",'.',"\t",txSt[key],"\t",txEnd[key],"\t","0,255,0\t",blocks[key],"\t",string.join(sz,','),"\t",string.join([str(i) for i in st],','), file=FO)
                                        else:
                                                print(key,"\t","pair-end mapped, but two ends mapped to different chromosome", file=FOF)
                                elif(count[key] >2):    #reads occur more than 2 times
                                        print(key,"\t","occurs more than 2 times in sam file", file=FOF)
                                        continue
                        FO.close()
                        FOF.close()
                        print("Done", file=sys.stderr)


        def getUnmap(self, outfile=None,fastq=True):
                '''Extract unmapped reads from SAM file and write to fastq [when fastq=True] or fasta [when fastq=False] file'''
                if outfile is None:
                        if fastq: outfile = self.fileName + ".unmap.fq"
                        else: outfile = self.fileName + ".unmap.fa"
                FO=open(outfile,'w')
                unmapCount=0
                print("Writing unmapped reads to\"",outfile,"\"... ", end=' ', file=sys.stderr)
                
                for line in self.f:
                        hits=[]
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip().split()     
                        flagCode=string.atoi(field[1])
                        seq=field[9]
                        qual=field[10]
                        if (flagCode & 0x0004) != 0:                                                    #read unmap
                                unmapCount +=1
                                if (flagCode & 0x0001) != 0:                                                    #paried in sequencing
                                        if (flagCode & 0x0040)!=0:seqID=field[0] + '/1'                 #first read     
                                        if (flagCode & 0x0080)!=0:seqID=field[0] + '/2'                 #second read
                                else: seqID=field[0]
                                
                                if fastq: FO.write('@' + seqID + '\n' + seq +'\n' + '+' +'\n' + qual+'\n')
                                else: FO.write('>' + seqID + '\n' + seq +'\n')

                print(str(unmapCount) + " reads saved!\n", file=sys.stderr)
                FO.close()
                self.f.seek(0)
        
        
        def getProperPair(self,outfile=None):
                '''Extract proper paried mapped reads.'''
                if outfile is None:
                        outfile = self.fileName + ".PP.sam"
                FO=open(outfile,'w')
                PPcount=0
                print("Writing proper paired reads to\"",outfile,"\"... ", end=' ', file=sys.stderr)
                for line in self.f:
                        hits=[]
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        if ((flagCode & 0x0001) != 0) and ((flagCode & 0x0002)!=0):
                                PPcount +=1
                                FO.write(line)
                FO.close()
                print(str(PPcount) + " reads were saved!\n", end=' ', file=sys.stderr)
                self.f.seek(0)

        def samNVC(self,outfile=None):
                '''for each read, calculate nucleotide frequency vs position'''
                if outfile is None:
                        outfile1 = self.fileName + ".NVC.xls"
                        outfile2 = self.fileName +".NVC_plot.r"
                else:
                        outfile1 = outfile + ".NVC.xls"
                        outfile2 = outfile +".NVC_plot.r"
                FO=open(outfile1,'w')
                RS=open(outfile2,'w')
                PPcount=0
                
                transtab = string.maketrans("ACGTNX","TGCANX")
                base_freq=collections.defaultdict(int)
                a_count=[]
                c_count=[]
                g_count=[]
                t_count=[]
                print("reading sam file ... ", file=sys.stderr)
                for line in self.f:
                        if line.startswith('@'):continue                                                                        #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue                                       #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        
                        if flagCode & 0x0010 ==0:       #plus strand
                                RNA_read = field[9].upper()
                        else:
                                RNA_read = field[9].upper().translate(transtab)[::-1]
                        for i in range(len(RNA_read)):
                                key = str(i) + RNA_read[i]
                                base_freq[key] += 1
                
                print("generating data matrix ...", file=sys.stderr)
                print("Position\tA\tC\tG\tT\tN\tX", file=FO)
                for i in range(len(RNA_read)):
                        print(str(i) + '\t', end=' ', file=FO)
                        print(str(base_freq[str(i) + "A"]) + '\t', end=' ', file=FO)
                        a_count.append(str(base_freq[str(i) + "A"]))
                        print(str(base_freq[str(i) + "C"]) + '\t', end=' ', file=FO)
                        c_count.append(str(base_freq[str(i) + "C"]))
                        print(str(base_freq[str(i) + "G"]) + '\t', end=' ', file=FO)
                        g_count.append(str(base_freq[str(i) + "G"]))
                        print(str(base_freq[str(i) + "T"]) + '\t', end=' ', file=FO)
                        t_count.append(str(base_freq[str(i) + "T"]))
                        print(str(base_freq[str(i) + "N"]) + '\t', end=' ', file=FO)
                        print(str(base_freq[str(i) + "X"]) + '\t', file=FO)
                FO.close()
                
                #generating R scripts
                print("generating R script  ...", file=sys.stderr)
                print("position=c(" + ','.join([str(i) for i in range(len(RNA_read))]) + ')', file=RS)
                print("A_count=c(" + ','.join(a_count) + ')', file=RS)
                print("C_count=c(" + ','.join(c_count) + ')', file=RS)
                print("G_count=c(" + ','.join(g_count) + ')', file=RS)
                print("T_count=c(" + ','.join(t_count) + ')', file=RS)
                print("total= A_count + C_count + G_count + T_count", file=RS)
                print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05", file=RS)
                print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total)", file=RS)
                
                print('pdf("NVC_plot.pdf")', file=RS)
                print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
                print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
                print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
                print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
                print('legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))', file=RS)
                print("dev.off()", file=RS)
                
                RS.close()
                #self.f.seek(0)

        def samGC(self,outfile=None):
                '''GC content distribution of reads'''
                if outfile is None:
                        outfile1 = self.fileName + ".GC.xls"
                        outfile2 = self.fileName +".GC_plot.r"
                else:
                        outfile1 = outfile + ".GC.xls"
                        outfile2 = outfile +".GC_plot.r"
                FO=open(outfile1,'w')
                RS=open(outfile2,'w')
                
                gc_hist=collections.defaultdict(int)    #key is GC percent, value is count of reads
                print("reading sam file ... ", file=sys.stderr)
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        gc_percent = "%4.2f" % ((field[9].upper().count('C') + field[9].upper().count('G'))/(len(field[9])+0.0)*100)
                        #print gc_percent
                        gc_hist[gc_percent] += 1
                
                print("writing GC content ...", file=sys.stderr)
                
                print("GC%\tread_count", file=FO)
                for i in list(gc_hist.keys()):
                        print(i + '\t' + str(gc_hist[i]), file=FO)
                        
                print("writing R script ...", file=sys.stderr)
                print("pdf('GC_content.pdf')", file=RS)
                print('gc=rep(c(' + ','.join([i for i in list(gc_hist.keys())]) + '),' + 'times=c(' + ','.join([str(i) for i in list(gc_hist.values())]) + '))', file=RS)
                print('hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100, file=RS)
                #print >>RS, "lines(density(gc),col='red')"
                print("dev.off()", file=RS)             
                #self.f.seek(0)
                
        def samDupRate(self,outfile=None,up_bound=500):
                '''Calculate reads's duplicate rates'''
                if outfile is None:
                        outfile1 = self.fileName + ".seq.DupRate.xls"
                        outfile2 = self.fileName + ".pos.DupRate.xls"
                        outfile3 = self.fileName +".DupRate_plot.r"
                else:
                        outfile1 = outfile + ".seq.DupRate.xls"
                        outfile2 = outfile + ".pos.DupRate.xls"
                        outfile3 = outfile +".DupRate_plot.r"
                SEQ=open(outfile1,'w')
                POS=open(outfile2,'w')
                RS=open(outfile3,'w')
                
                seqDup=collections.defaultdict(int)
                posDup=collections.defaultdict(int)
                
                seqDup_count=collections.defaultdict(int)
                posDup_count=collections.defaultdict(int)
                print("reading sam file ... ", file=sys.stderr)
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        if (flagCode & 0x0004) == 1: 
                                continue                                                                        #skip unmapped reads
                        seqDup[field[9]] +=1                    #key is read sequence
                        
                        #calculte duplicate read based on coordinates
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]       #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        chrom = field[2]
                        chromStart = string.atoi(field[3])-1
                        chromEnd=chromStart + sum(map(int,comb))
                        blockSize = []
                        for i in range(0,len(comb),2):
                                blockSize.append(str(comb[i]))
                        blockSizes = ','.join(blockSize)
                        blockStart=[]
                        for i in range(0,len(comb),2):
                                blockStart.append(str(sum(comb[:i])))
                        blockStarts = ','.join(blockStart)
                        
                        coord = chrom + ":" + str(chromStart) + "-" + str(chromEnd) + ":" + blockSizes + ":" + blockStarts
                        posDup[coord] +=1
                        
                print("report duplicte rate based on sequence ...", file=sys.stderr)
                print("Occurrence\tUniqReadNumber", file=SEQ)
                for i in list(seqDup.values()):                 #key is occurence, value is uniq reads number (based on seq)
                        seqDup_count[i] +=1
                for k in sorted(seqDup_count.keys()):   
                        print(str(k) +'\t'+ str(seqDup_count[k]), file=SEQ)
                SEQ.close()
                
                print("report duplicte rate based on mapping  ...", file=sys.stderr)
                print("Occurrence\tUniqReadNumber", file=POS)
                for i in list(posDup.values()):                 #key is occurence, value is uniq reads number (based on coord)
                        posDup_count[i] +=1
                for k in sorted(posDup_count.keys()):   
                        print(str(k) +'\t'+ str(posDup_count[k]), file=POS)
                POS.close()
                
                
                print("generate R script ...", file=sys.stderr)
                print("pdf('duplicateRead.pdf')", file=RS)
                print("par(mar=c(5,4,4,5),las=0)", file=RS)
                print("seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.keys()) ]) + ')', file=RS)
                print("seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.keys()) ]) + ')', file=RS)
                print("pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.keys()) ]) + ')', file=RS)
                print("pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.keys()) ]) + ')', file=RS)
                print("plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Frequency',pch=4,cex=0.8,col='blue',xlim=c(1,%d),yaxt='n')" % up_bound, file=RS)
                print("points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')", file=RS)
                print('ym=floor(max(log10(pos_uniqRead)))', file=RS)
                print("legend(%d,ym,legend=c('Sequence-base','Mapping-base'),col=c('red','blue'),pch=c(4,20))" % max(up_bound-200,1), file=RS)
                print('axis(side=2,at=0:ym,labels=0:ym)', file=RS)
                print('axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead)),round(pos_uniqRead[2]*100/sum(pos_uniqRead)),round(pos_uniqRead[3]*100/sum(pos_uniqRead)),round(pos_uniqRead[4]*100/sum(pos_uniqRead))))', file=RS)
                print('mtext(4, text = "Reads %", line = 2)', file=RS)
                #self.f.seek(0)
                
        def getUniqMapRead(self,outfile=None):
                '''Extract uniquely mapped reads.'''
                if outfile is None:
                        outfile = self.fileName + ".uniq.sam"
                FO=open(outfile,'w')
                Uniqcount=0
                print("Writing uniquely mapped reads to\"",outfile,"\"... ", end=' ', file=sys.stderr)
                for line in self.f:
                        hits=[]
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        if (flagCode & 0x0004) == 1: 
                                continue                        #skip unmapped reads
                        #else:
                                #print >>sys.stderr,line,
                        if (ParseSAM._uniqueHit_pat.search(line)):
                                print(line, end=' ', file=sys.stderr)
                                Uniqcount +=1
                                FO.write(line)
                FO.close()
                print(str(Uniqcount) + " reads were saved!\n", end=' ', file=sys.stderr)
                self.f.seek(0)

        def getWrongStrand(self,outfile=None):
                '''Extract pair-end reads mapped in incorrectly strand, such +/+ or -/-'''
                if outfile is None:
                        outfile = self.fileName + ".wrongStrand.sam"
                FO=open(outfile,'w')
                wrongStrand=0
                print("Writing incorrectly stranded reads to\"",outfile,"\"... ", end=' ', file=sys.stderr)
                for line in self.f:
                        hits=[]
                        if line.startswith('@'):continue                                                #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        if (flagCode & 0x0004) != 0: 
                                continue                        #skipped if read itself is unmapped
                        if (flagCode & 0x0008) !=0: 
                                continue                        #skipped if mate  is unmapped
                        if (flagCode & 0x0001) == 0: 
                                continue                        #skip single end sequencing'
                        if  ((flagCode & 0x0010) ==0) and ((flagCode & 0x0020)==0 ):
                                FO.write(line)
                                wrongStrand+=1
                        if  ((flagCode & 0x0010) !=0) and ((flagCode & 0x0020)!=0 ):
                                FO.write(line)
                                wrongStrand+=1

                FO.close()
                print(str(wrongStrand) + " reads were saved!\n", end=' ', file=sys.stderr)
                self.f.seek(0)
                
        def filterSpliceRead(self,outfile=None,min_overhang=8,min_gap=50,max_gap=1000000):
                '''filter spiced mapped reads from sam file. min_overhang is used to determine the reliability of splice sites
                splice reads with overhang size <8 will also be reported if the same splice sites has been suported by
                at least 1 read with overhang size >8. Multiple spliced reads (belong to the same splice junction) 
                will always be reported. min_overhang, min_gap and max_gap only applied to one-time splice read'''
                
                if outfile is None:
                        outfile = self.fileName + ".SR.sam"
                        #outfile2 = self.fileName + ".SR.filter.sam"                    
                splice_sites=collections.defaultdict(set)
                print("\tDetermine splice sites with proper overhang, intron size ... ", end=' ', file=sys.stderr)
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        if not (ParseSAM._uniqueHit_pat.search(line)):  #skip non unique mapped read
                                continue
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        map_st = int(field[3])
                        chrom = field[2]
                        if (flagCode & 0x0004) == 1:continue                    #skip unmapped reads
                        
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]       #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        if (len(comb)==1):      #skip non-spliced
                                #print line,
                                continue
                        if (len(comb)>3):       #skip multiple spliced
                                continue
                        else:                           #one-time spliced
                                if (comb[1] < min_gap or comb[1] > max_gap):
                                        continue
                                else:
                                        if (comb[0] >= min_overhang):
                                                splice_sites[chrom].add(map_st + comb[0])
                                        if (comb[2] >= min_overhang):
                                                splice_sites[chrom].add(map_st + comb[0] + comb[1])                     
                self.f.seek(0)  
                print("Done", file=sys.stderr)
                
                
                FO=open(outfile,'w')
                #FO2=open(outfile2,'w')
                print("\tExtracting splicing reads  ... ", end=' ', file=sys.stderr)
                total_SR =0
                extract_SR =0
                total_read =0
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        map_st = int(field[3])
                        chrom = field[2]
                        if (flagCode & 0x0004) == 1:continue                    #skip unmapped reads
                        total_read +=1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]       #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        if (len(comb)==1):      #skip non-spliced
                                continue
                        total_SR +=1
                        
                        
                        if (len(comb)>3):       #multipel splice read. report directly
                                FO.write(line)
                                extract_SR +=1
                        else:                           #one-time splice read
                                if (comb[1] < min_gap or comb[1] > max_gap):
                                        continue
                                if (chrom in splice_sites) and ((map_st + comb[0]) in splice_sites[chrom]) and ((map_st + comb[0] + comb[1]) in splice_sites[chrom]):
                                        FO.write(line)
                                        #print line
                                        extract_SR +=1
                                else:
                                        #FO2.write(line)
                                        continue
                print("Done", file=sys.stderr)
                print("\tTotal mapped Read: " + str(total_read), file=sys.stderr)
                print("\tTotal Splicing Read: " + str(total_SR), file=sys.stderr)
                print("\\Usable Splicing Read: " + str(extract_SR), file=sys.stderr)
                FO.close()
                #FO2.close()
                self.f.seek(0)  

        def getSpliceRead(self,outfile=None):
                '''Extract spiced mapped reads from sam file'''
                
                if outfile is None:
                        outfile = self.fileName + ".SR.sam"
                FO=open(outfile,'w')
                print("\tExtract splicing reads without any filter ...", end=' ', file=sys.stderr)
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines
                        field=line.rstrip('\n').split() 
                        flagCode=string.atoi(field[1])
                        if (flagCode & 0x0004) == 1:continue                    #skip unmapped reads                    
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]       #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        if (len(comb)>=3):
                                FO.write(line)
                
                print("Done", file=sys.stderr)
                self.f.seek(0)  
                FO.close()

        def collapseSAM(self, outfile=None,collevel=10):
                '''At most collevel[default=10] identical reads will be retained in outputting SAM file
                The original SAM file must be sorted before hand. if not, using linux command like "sort -k3,3 -k4,4n myfile.sam >myfile.sorted.sam" '''
                if outfile is None:
                        outfile = self.fileName + ".collapsed.sam"
                print("Writing collapsed SAM file to\"",outfile,"\"... ", file=sys.stderr)
                FO=open(outfile,'w')            
                flag=""
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines       
                        field=line.rstrip('\n').split() 
                        if (string.atoi(field[1]) & 0x0004)!=0: continue        #skip unmapped line     
                        id=field[2] + field[3] + field[5]
                        if (id != flag):
                                FO.write(line)
                                flag=id
                                skipTrigger=0
                        else:
                                skipTrigger+=1
                                if skipTrigger < collevel:
                                        FO.write(line)
                                else:continue
                FO.close()
                self.f.seek(0)

        def qualSAM(self,read_len,outfile=None):
                '''calculate phred quality score for each base in read (5->3)'''
                if outfile is None:
                        outfile = self.fileName + ".qual.plot.r"
                else:
                        outfile = outfile + ".qual.plot.r"
                FO=open(outfile,'w')
                print("\tcalculating quality score ... ", file=sys.stderr)
                qual_min={}
                qual_max={}
                qual_sum={}
                total_read=0
                for i in range(0,read_len):
                        qual_min.setdefault(i,1000)
                        qual_max.setdefault(i,-1)
                        qual_sum.setdefault(i,0.0)
                
                for line in self.f:
                        if line[0] == '@':continue                                              #skip head lines        
                        if ParseSAM._reExpr2.match(line):continue               #skip blank lines       
                        field=line.rstrip('\n').split()
                        #if (string.atoi(field[1]) & 0x0004)!=0: continue       #skip unmapped line
                        
                        if (len(field[10]) != read_len):
                                continue
                        if (string.atoi(field[1]) & 0x0010)==0: #query map to +
                                qual_str=field[10]
                        else:
                                qual_str=field[10][::-1]
                        total_read +=1
                        for i in range(0,read_len):
                                #print ord(qual_str[i])-33,
                                qual_sum[i] += ord(qual_str[i])-33
                                if(qual_min[i] > (ord(qual_str[i])-33)):
                                        qual_min[i] = ord(qual_str[i])-33
                                if(qual_max[i] < (ord(qual_str[i])-33)):
                                        qual_max[i] = ord(qual_str[i])-33
                        #print '\n',
                min_qualities = [str(qual_min[i]) for i in range(0,read_len)]
                max_qualities =[str(qual_max[i]) for i in range(0,read_len)]
                avg_qualities = [str(qual_sum[i]/total_read) for i in range(0,read_len)]
                nt_pos = [str(i) for i in range(0,read_len)]
                print("nt_pos=c(" + ','.join(nt_pos)  + ')', file=FO)
                print("max_qual=c(" + ','.join(max_qualities)  + ')', file=FO)
                print("min_qual=c(" + ','.join(min_qualities)  + ')', file=FO)
                print("avg_qual=c(" + ','.join(avg_qualities)  + ')', file=FO)
                print("pdf('phred_qual.pdf')", file=FO)
                print("plot(nt_pos,avg_qual, xlab=\"Nucleotide Position (5'->3')\", ylab='Phred Quality',ylim=c(0,97),lwd=2,type='s')", file=FO)
                print('lines(nt_pos,max_qual,type="s",lwd=2,col="red")', file=FO)
                print('lines(nt_pos,min_qual,type="s",lwd=2,col="blue")', file=FO)
                print('legend(0,100,legend=c("Max","Average","Min"),col=c("red","black","blue"),lwd=2)', file=FO)
                print('dev.off()', file=FO)
                #for i in range(0,read_len):
                #       print >>sys.stderr, str(i) + '\t' + str(qual_max[i]) + '\t' + str(qual_min[i]) + '\t' + str(qual_sum[i]/total_read)
                #self.f.seek(0)


        def samToBinnedArray(self):
                """Convert SAM file to BinnedArray."""
                
                lines=0
                for line in self.f:
                        if line.startswith('@'):continue                                                #skip head lines
                        if ParseSAM._reExpr2.match(line):continue                               #skip blank lines
                        field=line.rstrip().split()     
                        if (string.atoi(field[1]) & 0x0004)!=0: continue                #skip unmapped line             
                        txStart=string.atoi(field[3])
                        #if (string.atoi(field[1]) & 0x0010 != 0):
                        #       strand='-'
                        #else:
                        #       strand='+'
                        lines +=1
                        scores={}
                        chrom = field[2]
                        
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(field[5])]       #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        if not chrom in scores:scores[chrom] = BinnedArray()
                        
                        for i in range(0,len(comb),2):
                                for pos in range(txStart + sum(comb[:i]), txStart + sum(comb[:i]) + comb[i]):
                                        tmp = scores[chrom][pos]
                                        if isNaN(tmp):                          
                                                scores[chrom][pos] =1
                                        else:
                                                scores[chrom][pos] +=1                  
                        if lines % 10000 == 0: print("%i lines loaded \r" % lines, file=sys.stderr)
                return scores
                self.f.seek(0)

class QCSAM(object):
        '''Perform basic quality control. Useful for RNA-seq experiment'''
        

        def __init__(self,samFile):
                '''constructor'''
                if samFile == '-':
                        self.fileName = "STDIN"
                        self.f = sys.stdin
                else:
                        self.fileName = os.path.basename(samFile)
                        self.f=open(samFile,'r')

                
        def distribSAM(self,refbed,outfile=None):
                '''calculate reads distribution over genome features (Exon reads, Inton reads, Intergenic 
                reads, spliced reads). A bed file representing the gene model (i.e. refseq) must be provided
                Two bed format files will be generated: outfile_exon.count.bed and outfile_intron.count.bed.
                The 5th column is number of reads fallen into the region defined by the first 3 columns'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if outfile is None:
                        exon_count = self.fileName + "_exon.count.bed"
                        intron_count = self.fileName + "_intron.count.bed"
                        rscript=self.fileName + ".piechart.r"
                        rpdf=self.fileName + ".piechart.pdf"
                else:
                        exon_count = outfile + "_exon.count.bed"
                        intron_count = outfile +  "_intron.count.bed"
                        rscript=outfile + ".piechart.r"
                        rpdf=outfile + ".piechart.pdf"
                
                EXON_OUT = open(exon_count,'w')
                INTRON_OUT =open(intron_count,'w')
                R_OUT = open(rscript,'w')
                
                ranges={}
                intronReads=0
                exonReads=0
                intergenicReads=0
                totalReads=0
                splicedReads=0
                
                #read SAM 
                print("reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue           #skip unmap reads
                        totalReads +=1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        if( len(comb)>1):
                                splicedReads +=1        
                                continue
                        else:
                                chrom=fields[2].upper()
                                #st=int(fields[3])-1
                                #end= st +len(fields[9])
                                mid = int(fields[3]) + int(len(fields[9])/2)
                                if chrom not in ranges:
                                        ranges[chrom] = Intersecter()
                                ranges[chrom].add_interval( Interval( mid, mid ) )
                        
                self.f.seek(0)
                print("Done", file=sys.stderr)
                
                #read refbed file
                print("Assign reads to "+ refbed + '...', end=' ', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith('#'):continue
                                if line.startswith('track'):continue
                                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                                intron_starts = exon_ends[:-1]
                                intron_ends=exon_starts[1:]
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue

                                # assign reads to intron                                
                        if(strand == '-'):
                                intronNum=len(intron_starts)
                                exonNum=len(exon_starts)
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                intronReads += hits
                                                INTRON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\n')
                                                intronNum -= 1
                                                        
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                exonReads += hits
                                                EXON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\n')
                                                exonNum -= 1
                        elif(strand == '+'):
                                intronNum=1
                                exonNum=1
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                intronReads += hits
                                                INTRON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\n')
                                                intronNum += 1    
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                exonReads += hits
                                                EXON_OUT.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\n')
                                                exonNum += 1            
                intergenicReads=totalReads-exonReads-intronReads-splicedReads           
                print("Done." + '\n', file=sys.stderr)
                print("Total reads:\t" + str(totalReads), file=sys.stderr)
                print("Exonic reads:\t" + str(exonReads), file=sys.stderr) 
                print("Intronic reads:\t" + str(intronReads), file=sys.stderr) 
                print("Splicing reads:\t" + str(splicedReads), file=sys.stderr)
                print("Intergenic reads:\t" + str(intergenicReads), file=sys.stderr)
                
                print("writing R script ...", end=' ', file=sys.stderr)
                totalReads=float(totalReads)
                print("pdf('%s')" % rpdf, file=R_OUT)
                print("dat=c(%d,%d,%d,%d)" % (exonReads,splicedReads,intronReads,intergenicReads), file=R_OUT)
                print("lb=c('exon(%.2f)','junction(%.2f)','intron(%.2f)','intergenic(%.2f)')" % (exonReads/totalReads,splicedReads/totalReads,intronReads/totalReads,intergenicReads/totalReads), file=R_OUT)
                print("pie(dat,labels=lb,col=rainbow(4),clockwise=TRUE,main='Total reads = %d')" % int(totalReads), file=R_OUT)
                print("dev.off()", file=R_OUT)
                print("Done.", file=sys.stderr)
        
        
        def coverageGeneBody(self,refbed,outfile=None):
                '''Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divied
                into 100 regsions'''
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if outfile is None:
                        outfile1 = self.fileName + ".geneBodyCoverage_plot.r"
                        outfile2 =  self.fileName + ".geneBodyCoverage.txt"
                else:
                        outfile1 = outfile + ".geneBodyCoverage_plot.r"
                        outfile2 = outfile + ".geneBodyCoverage.txt"
                OUT1 = open(outfile1,'w')
                OUT2 = open(outfile2,'w')

                ranges={}
                totalReads=0
                fragment_num=0          #splice reads will counted twice
                rpkm={}
                
                #read SAM 
                print("reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue           #skip unmap reads
                        totalReads +=1
                        
                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        fragment_num += (len(comb) +1)/2
                        blockStart=[]
                        blockSize=[]
                        
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        for st,size in zip(blockStart,blockSize):
                                if chrom not in ranges:
                                        ranges[chrom] = Intersecter()
                                ranges[chrom].add_interval( Interval( st, st+size ) )
                print("Done", file=sys.stderr)          

                print("calculating coverage over gene body ...", file=sys.stderr)
                coverage=collections.defaultdict(int)
                flag=0
                for line in open(refbed,'r'):
                        try:
                                if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5]
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue
                        gene_all_base=[]
                        percentile_base=[]
                        mRNA_len =0
                        flag=0
                        for st,end in zip(exon_starts,exon_ends):
                                gene_all_base.extend(list(range(st+1,end+1)))           #0-based coordinates on genome
                                mRNA_len = len(gene_all_base)
                                if mRNA_len <100:
                                        flag=1
                                        break
                        if flag==1: continue
                        if strand == '-':
                                gene_all_base.sort(reverse=True)                        #deal with gene on minus stand
                        else:
                                gene_all_base.sort(reverse=False)
                        percentile_base = mystat.percentile_list (gene_all_base)        #get 101 points from each gene's coordinates
                        
                        for i in range(0,len(percentile_base)):
                                if chrom in ranges:
                                        coverage[i] += len(ranges[chrom].find(percentile_base[i], percentile_base[i]+1))
                x_coord=[]
                y_coord=[]
                print("Total reads: " + str(totalReads), file=OUT2)
                print("Fragment number: " + str(fragment_num), file=OUT2)
                print("percentile\tcount", file=OUT2)
                for i in coverage:
                        x_coord.append(str(i))
                        y_coord.append(str(coverage[i]))
                        print(str(i) + '\t' + str(coverage[i]), file=OUT2)
                print("pdf('geneBody_coverage.pdf')", file=OUT1)
                print("x=0:100", file=OUT1)
                print("y=c(" + ','.join(y_coord) + ')', file=OUT1)
                print("plot(x,y,xlab=\"percentile of gene body (5'->3')\",ylab='read number',type='s')", file=OUT1)
                print("dev.off()", file=OUT1)
                        
        def calculateRPKM(self,refbed,outfile=None):
                '''calculate RPKM values for each gene in refbed. Only uniquely aligned reads are used. 
                Spilced reads are split. output raw read connt and eRPKM (eRPKM = exon Represented times Per Kb 
                exon per Million mapped reads) for each exon, intron and mRNA'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if outfile is None:
                        rpkm_file = self.fileName + ".rpkm.xls"
                else:
                        rpkm_file = outfile + ".rpkm.xls"
                RPKM_OUT=open(rpkm_file,'w')
                
                ranges={}
                totalReads=0
                cUR=0
                sR=0
                multiMapReads=0
                rpkm={}
                
                #read SAM 
                print("reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue           #skip unmap reads
                        totalReads +=1
                        if not ParseSAM._uniqueHit_pat.search(line):            #skip multiple mapped reads
                                multiMapReads +=1
                                continue

                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        cUR += (len(comb) +1)/2
                        if(len(comb)>1):
                                sR+=1
                        blockStart=[]
                        blockSize=[]
                        
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        for st,size in zip(blockStart,blockSize):
                                mid = int(st) + (size/2)
                                if chrom not in ranges:
                                        ranges[chrom] = Intersecter()
                                ranges[chrom].add_interval( Interval( mid, mid ) )
                        
                self.f.seek(0)
                print("Done", file=sys.stderr)
                print("Total mapped reads (TR): " + str(totalReads), file=RPKM_OUT) 
                print("Multiple mapped reads (MR): " + str(multiMapReads), file=RPKM_OUT)
                print("Uniquely mapped reads (UR): " + str(totalReads - multiMapReads), file=RPKM_OUT)
                print("Spliced  mapped reads (SR): " + str(sR), file=RPKM_OUT)
                print("Corrected uniquely mapped reads (cUR): " + str(cUR), file=RPKM_OUT)
                if totalReads ==0:
                        sys.exit(1)
                
                #read refbed file
                print("Assign reads to "+ refbed + '...', end=' ', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith('#'):continue
                                if line.startswith('track'):continue
                                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                intron_starts = exon_ends[:-1]
                                intron_ends=exon_starts[1:]
                                key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue

                        # assign reads to intron                                
                        mRNA_count=0
                        mRNA_len=sum(exon_sizes)
                        if(strand == '-'):
                                intronNum=len(intron_starts)
                                exonNum=len(exon_starts)
                                                                
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                intronNum -= 1
                                                        
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))                                   
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                exonNum -= 1
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
                                        rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                                        rpkm[key] = 0
                        elif(strand == '+'):
                                intronNum=1
                                exonNum=1
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                intronNum += 1    
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                exonNum += 1            
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
                                        rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                                        rpkm[key] = 0
                print("Done", file=sys.stderr)
                return rpkm
                self.f.seek(0)

        def calculateRPKM2(self,refbed,outfile=None):
                '''calculate RPKM values for each gene in refbed. Only uniquely aligned reads are used. 
                Spilced reads are split. output raw read connt and eRPKM (eRPKM = exon Represented times Per Kb 
                exon per Million mapped reads) for each exon, intron and mRNA
                NOTE: intronic reads are not counted as part of total reads'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if outfile is None:
                        rpkm_file = self.fileName + ".rpkm.xls"
                else:
                        rpkm_file = outfile + ".rpkm.xls"
                RPKM_OUT=open(rpkm_file,'w')
                
                ranges={}
                exon_ranges={}
                totalReads=0
                #intronic=0
                cUR=0
                sR=0
                multiMapReads=0
                rpkm={}
                
                #read gene model file, the purpose is to remove intronic reads
                print("Reading reference gene model "+ refbed + '...', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith(('#','track','browser')):continue
 
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue        

                        for st,end in zip(exon_starts,exon_ends):                               
                                if chrom not in exon_ranges:
                                        exon_ranges[chrom] = Intersecter()
                                exon_ranges[chrom].add_interval( Interval( st, end ) )          

                #read SAM 
                print("reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue           #skip unmap reads
                        totalReads +=1
                        if not ParseSAM._uniqueHit_pat.search(line):            #skip multiple mapped reads
                                multiMapReads +=1
                                continue

                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        #cUR += (len(comb) +1)/2
                        if(len(comb)>1):
                                sR+=1
                        blockStart=[]
                        blockSize=[]
                        
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        #build bitset only for exonic reads
                        for st,size in zip(blockStart,blockSize):
                                if (chrom in exon_ranges) and (len(exon_ranges[chrom].find(st,st+size)) >0):    #if we found this fragment is overlapped with exon
                                        cUR += 1
                                        mid = int(st) + (size/2)
                                        if chrom not in ranges:
                                                ranges[chrom] = Intersecter()
                                        ranges[chrom].add_interval( Interval( mid, mid ) )
                        else:                                                                                                                                                   #if this framgnet is intronic, skip it.
                                #intronic +=1
                                continue        
                self.f.seek(0)
                print("Done", file=sys.stderr)
                print("Total mapped reads (TR): " + str(totalReads), file=RPKM_OUT) 
                print("Multiple mapped reads (MR): " + str(multiMapReads), file=RPKM_OUT)
                print("Uniquely mapped reads (UR): " + str(totalReads - multiMapReads), file=RPKM_OUT)
                print("Spliced  mapped reads (SR): " + str(sR), file=RPKM_OUT)
                print("Corrected uniquely mapped reads (cUR, non-intronic fragments): " + str(cUR), file=RPKM_OUT)
                #print >>RPKM_OUT, "Intronic Fragments (IF): " + str(intronic)
                if totalReads ==0:
                        sys.exit(1)
                
                #read refbed file
                print("Assign reads to "+ refbed + '...', end=' ', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith('#'):continue
                                if line.startswith('track'):continue
                                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                intron_starts = exon_ends[:-1]
                                intron_ends=exon_starts[1:]
                                key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue

                        # assign reads to intron                                
                        mRNA_count=0
                        mRNA_len=sum(exon_sizes)
                        if(strand == '-'):
                                intronNum=len(intron_starts)
                                exonNum=len(exon_starts)
                                                                
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                intronNum -= 1
                                                        
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))                                   
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                exonNum -= 1
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
                                        rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                                        rpkm[key] = 0
                        elif(strand == '+'):
                                intronNum=1
                                exonNum=1
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                intronNum += 1    
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                hits= len(ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(cUR))) +'\n')
                                                exonNum += 1            
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*cUR)) +'\n')
                                        rpkm[key] = mRNA_count*1000000000.0/(mRNA_len*cUR)
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                                        rpkm[key] = 0
                print("Done", file=sys.stderr)
                return rpkm
                self.f.seek(0)

        def filterKnownReads(self,refbed,outfile=None):
                '''Compare SAM files with reference gene model, all reads mapped to gene model will be filted
                out. The remainning unknown reads will be writtern to a new SAM file'''
                
                totalReads=0    #total mapped reads
                unknownReads=0
                ranges={}
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                
                if outfile is None:
                        out_file = self.fileName + ".unknownReads.SAM"
                else:
                        out_file = outfile + ".unknownReads.SAM"        
                OUT=open(out_file,'w')
                
                print("Reading reference gene model "+ refbed + '...', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith(('#','track','browser')):continue
 
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue        

                        for st,end in zip(exon_starts,exon_ends):                               
                                if chrom not in ranges:
                                        ranges[chrom] = Intersecter()
                                ranges[chrom].add_interval( Interval( st, end ) )               

                print("Processing SAM file "+ self.fileName + '...', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue                   #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):    #skip multiple mapped reads
                                continue
                        
                        blockStart=[]
                        blockSize=[]
                        totalReads +=1
                        
                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']                    
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])

                        for st,size in zip(blockStart,blockSize):
                                if (chrom in ranges) and (len(ranges[chrom].find(st,st+size)) >0):      #if we found this read is overlapped with known gene
                                        break
                        else:
                                OUT.write(line)
                                unknownReads +=1
                OUT.close()
                print("Total reads mapped to genome: " + str(totalReads), file=sys.stderr)
                print("Total reads not overlapped with any exon: " + str(unknownReads), file=sys.stderr)
                self.f.seek(0)

        def genomicFragSize(self,outfile=None,low_bound=0,up_bound=1000,step=10):
                '''estimate the genomic fragment size of mRNA experiment. fragment size = insert_size + 2 x read_length'''
                
                
                if outfile is None:
                        out_file1 = self.fileName + ".fragSize.txt"
                        out_file2 = self.fileName + ".fragSize.Freq.txt"
                        out_file3 = self.fileName + ".fragSize_plot.r"
                else:
                        out_file1 = outfile + ".fragSize.txt"   
                        out_file2 = outfile + ".fragSize.Freq.txt"
                        out_file3 = outfile + ".fragSize_plot.r"
                
                FO=open(out_file1,'w')
                FQ=open(out_file2,'w')
                RS=open(out_file3,'w')
                
                chrom="chr100"          #this is the fake chromosome
                ranges={}
                ranges[chrom]=Intersecter()
                
                window_left_bound = list(range(low_bound,up_bound,step))
                frag_size=0

                pair_num=0.0
                ultra_low=0.0
                ultra_high=0.0
                size=[]
                counts=[]
                count=0
                print("Reading SAM file "+ self.fileName + ' ... ', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        #if fields[0] in pairRead_info:
                        #       continue
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0001) ==0:
                                print("NOT pair-end sequencing", file=sys.stderr)
                                sys.exit(1)
                        if (flagCode & 0x0004) != 0: continue                   #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):    #skip multiple mapped reads
                                continue
                        if (flagCode & 0x0008 !=0):                                             #skip single-end mapped reads
                                continue
                        if (fields[7] =='0'):
                                continue
                        if (int(fields[3]) > int(fields[7])):                   #left < right
                                continue
                        pair_num +=1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        read_len = len(fields[9])
                        if (len(comb)==1):              # this read is NOT spliced
                                frag_size = (int(fields[7]) - int(fields[3]) +1) + read_len
                        elif (len(comb) >1):    # this read is spliced
                                frag_size = (int(fields[7]) - int(fields[3]) +1) + read_len - sum(comb[1::2])
                        FO.write(fields[0] + '\t' + str(frag_size) + '\n')
                        if frag_size <= low_bound:
                                ultra_low+=1
                                continue
                        elif frag_size > up_bound:
                                ultra_high +=1
                                continue
                        ranges[chrom].add_interval( Interval( frag_size-1, frag_size ) )
                print("Done", file=sys.stderr)
                
                if pair_num==0:
                        print("Cannot find paired reads", file=sys.stderr)
                        sys.exit(0)
                print("Total paired read " + str(pair_num), file=FQ)
                print("<=" + str(low_bound) + "\t"+ str(ultra_low), file=FQ)
                for st in window_left_bound:
                        size.append(str(st + step/2))
                        count = str(len(ranges[chrom].find(st,st + step)))
                        counts.append(count)
                        print(str(st) + '\t' + str(st+step) +'\t' + count, file=FQ)             
                print(">" + str(up_bound) + "\t"+ str(ultra_high), file=FQ)
                
                print("pdf('gFragSize.pdf')", file=RS)
                print("par(mfrow=c(2,1),cex.main=0.8,cex.lab=0.8,cex.axis=0.8,mar=c(4,4,4,1))", file=RS)
                print('pie(c(%d,%d,%d),col=rainbow(3),cex=0.5,radius=1,main="Total %d fragments",labels=c("fraSize <= %d\\n(%4.2f%%)","fragSize > %d\\n(%4.2f%%)","%d < fragSize <= %d\\n(%4.2f%%)"), density=rep(80,80,80),angle=c(90,140,170))' % (ultra_low, ultra_high, pair_num -ultra_low -ultra_high, pair_num, low_bound, ultra_low*100/pair_num, up_bound, ultra_high*100/pair_num, low_bound, up_bound, 100-ultra_low*100/pair_num - ultra_high*100/pair_num), file=RS)
                print('fragsize=rep(c(' + ','.join(size) + '),' + 'times=c(' + ','.join(counts) + '))', file=RS)
                print('frag_sd = round(sd(fragsize))', file=RS)
                print('frag_mean = round(mean(fragsize))', file=RS)
                print('hist(fragsize,probability=T,breaks=%d,xlab="Fragment size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(window_left_bound), file=RS)
                print("lines(density(fragsize,bw=%d),col='red')" % (2*step), file=RS)
                print("dev.off()", file=RS)
                FO.close()
                FQ.close()
                RS.close()
                #self.f.seek(0)
                
                
        def saturation_RPKM(self,refbed,outfile=None,sample_start=5,sample_step=5,sample_end=100):
                '''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if outfile is None:
                        rpkm_file = self.fileName + ".eRPKM.xls"
                        raw_file = self.fileName + ".rawCount.xls"
                else:
                        rpkm_file = outfile + ".eRPKM.xls"
                        raw_file = outfile + ".rawCount.xls"
                
                RPKM_OUT = open(rpkm_file,'w')
                RAW_OUT = open(raw_file ,'w')
                
                ranges={}
                totalReads=0
                cUR_num = 0     #number
                block_list=[]   #non-spliced read AS IS, splicing reads were counted multiple times
                                
                #read SAM 
                my_pat = re.compile(r'NH:i:(\d+)\b')
                NH_tag=0
                print("Reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue           #skip unmap reads
                        totalReads +=1
                        hitNum =[int(i) for i in my_pat.findall(line)]
                        if len(hitNum) ==0:
                                NH_tag=1                                                                #cannot determine uniqness without NH tag
                        elif len(hitNum) ==1:
                                if int(hitNum[0])>1: continue                           #skip multiple mapped reads
                        else:
                                print("More than 1 NH tag found within a single line. Incorrect SAM format!", file=sys.stderr)
                                sys.exit(1)

                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        cUR_num += (len(comb) +1)/2
                        blockStart=[]
                        blockSize=[]
                        
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        for st,size in zip(blockStart,blockSize):
                                mid = int(st) + (size/2)
                                block_list.append(chrom + ":" + str(mid))
                
                if NH_tag==1:
                        print("Warn: NO NH tag found. Cannot determine uniqueness of alignment. All alignments will be used", file=sys.stderr)
                print("Done", file=sys.stderr)
                
                print("shuffling alignments ...", end=' ', file=sys.stderr)
                random.shuffle(block_list)
                print("Done", file=sys.stderr)
                
                
                ranges={}
                sample_size=0
                frag_total = cUR_num
                RPKM_table=collections.defaultdict(list)
                rawCount_table=collections.defaultdict(list)
                RPKM_head=['chr','start','end','name','score','strand']

                tmp=list(range(sample_start,sample_end,sample_step))
                tmp.append(100)
                #=========================sampling uniquely mapped reads from population
                for pertl in tmp:       #[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
                        index_st = int(frag_total * (pertl-sample_step)/100.0)
                        index_end = int(frag_total * pertl/100.0)
                        if index_st < 0: index_st = 0
                        sample_size += index_end -index_st
                        
                        RPKM_head.append(str(pertl) + '%')
                        print("sampling " + str(pertl) +"% (" + str(sample_size) + ") fragments ...", end=' ', file=sys.stderr)
                        for i in range(index_st, index_end):
                                (chr,coord) = block_list[i].split(':')
                                if chr not in ranges:
                                        ranges[chr] = Intersecter()
                                ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                
                        #========================= calculating RPKM based on sub-population
                        #print >>sys.stderr, "assign reads to "+ refbed + '...',
                        for line in open(refbed,'r'):
                                try:
                                        if line.startswith(('#','track','browser')):continue  
                        # Parse fields from gene tabls
                                        fields = line.split()
                                        chrom     = fields[0].upper()
                                        tx_start  = int( fields[1] )
                                        tx_end    = int( fields[2] )
                                        geneName      = fields[3]
                                        strand    = fields[5].replace(" ","_")
                                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                        exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                        key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                                except:
                                        print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                                        continue
                                mRNA_count=0    #we need to initializ it to 0 for each gene
                                mRNA_len=sum(exon_sizes)
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in ranges:
                                                mRNA_count += len(ranges[chrom].find(st,end))           
                                if mRNA_len ==0:
                                        print(geneName + " has 0 nucleotides. Exit!", file=sys.stderr)
                                        sys.exit(1)
                                if sample_size == 0:
                                        print("Too few reads to sample. Exit!", file=sys.stderr)
                                        sys.exit(1)
                                mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * sample_size)
                                RPKM_table[key].append(str(mRNA_RPKM))
                                rawCount_table[key].append(str(mRNA_count))
                        print("Done", file=sys.stderr)

                #self.f.seek(0)
                print('\t'.join(RPKM_head), file=RPKM_OUT)
                print('\t'.join(RPKM_head), file=RAW_OUT)
                for key in RPKM_table:
                        print(key + '\t', end=' ', file=RPKM_OUT)
                        print('\t'.join(RPKM_table[key]), file=RPKM_OUT)
                        print(key + '\t', end=' ', file=RAW_OUT)
                        print('\t'.join(rawCount_table[key]), file=RAW_OUT)                     

        def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1):
                '''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction'''
                
                if outfile is None:
                        out_file = self.fileName + ".junctionSaturation_plot.r"
                else:
                        out_file = outfile + ".junctionSaturation_plot.r"
                if refgene is None:
                        print("You must provide reference gene model in bed format.", file=sys.stderr)
                        sys.exit(1)
                
                OUT = open(out_file,'w')


                #reading reference gene 
                knownSpliceSites= set()
                print("reading reference bed file: ",refgene, " ... ", end=' ', file=sys.stderr)
                for line in open(refgene,'r'):
                        if line.startswith(('#','track','browser')):continue  
                        fields = line.split()
                        if(len(fields)<12):
                                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                                continue
                        chrom     = fields[0].upper()
                        tx_start = int( fields[1] )
                        tx_end   = int( fields[2] )
                        if int(fields[9] ==1):
                                continue        
                        
                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        intron_start = exon_ends[:-1]
                        intron_end=exon_starts[1:]
                        for st,end in zip (intron_start, intron_end):
                                knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))
                print("Done! Total "+str(len(knownSpliceSites)) + " known splicing sites", file=sys.stderr)


                #read SAM file
                samSpliceSites=[]
                intron_start=[]
                intron_end=[]
                uniqSpliceSites=collections.defaultdict(int)
                print("Reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        
                        if (flagCode & 0x0004) != 0: continue                           #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):            #skip multiple mapped reads
                                continue
                        if (len(ParseSAM._splicedHit_pat.findall(fields[5]))==1):       #skip non-spilced reads
                                continue

                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        blockStart=[]
                        blockSize=[]
                        blockEnd=[]             
                        #if intron size < min_intron, skip. 
                        flag=0
                        for i in range(1,len(comb),2):
                                if comb[i] < min_intron:
                                        flag=1
                                        break
                        if flag ==1:
                                continue                
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        for st,size in zip(blockStart,blockSize):
                                end = st + size
                                blockEnd.append(end)
                        intron_st = blockEnd[:-1]
                        intron_end = blockStart[1:]
                        
                        for st,end in zip(intron_st, intron_end):
                                samSpliceSites.append(chrom + ":" + str(st) + "-" + str(end))                   
                #self.f.seek(0)
                print("Done", file=sys.stderr)
                


                print("shuffling alignments ...", end=' ', file=sys.stderr)
                random.shuffle(samSpliceSites)
                print("Done", file=sys.stderr)
                                
                #resampling
                SR_num = len(samSpliceSites)
                sample_size=0
                knownSpliceSites_num = 0
                known_junc=[]
                all_junc=[]
                #=========================sampling uniquely mapped reads from population
                tmp=list(range(sample_start,sample_end,sample_step))
                tmp.append(100)
                for pertl in tmp:       #[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
                        knownSpliceSites_num = 0
                        index_st = int(SR_num * ((pertl - sample_step)/100.0))
                        index_end = int(SR_num * (pertl/100.0))
                        if index_st < 0: index_st = 0
                        sample_size += index_end -index_st
                        
                        print("sampling " + str(pertl) +"% (" + str(sample_size) + ") unique splicing alignments ...", end=' ', file=sys.stderr)
                        for i in range(index_st, index_end):
                                uniqSpliceSites[samSpliceSites[i]] +=1  
                        all_junc.append(str(len(list(uniqSpliceSites.keys()))))
                        for sj in uniqSpliceSites:
                                if sj in knownSpliceSites and uniqSpliceSites[sj] >= recur:
                                        knownSpliceSites_num +=1
                        print(str(knownSpliceSites_num) + " known splicing junctions", file=sys.stderr)
                        known_junc.append(str(knownSpliceSites_num))
                        
                #for j in uniq_SJ:
                        #print >>OUT, j + "\t" + str(uniq_SJ[j])
                print("pdf('junction_saturation.pdf')", file=OUT)
                print("x=c(" + ','.join([str(i) for i in tmp]) + ')', file=OUT)
                print("y=c(" + ','.join(known_junc) + ')', file=OUT)
                print("z=c(" + ','.join(all_junc) + ')', file=OUT)
                print("plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(%d,%d))" % (int(int(known_junc[0])/1000), int(int(all_junc[-1])/1000)), file=OUT)
                print("points(x,y/1000,type='o',col='red')", file=OUT)
                print('legend(5,%d, legend=c("All detected junction","Annotated junction"),col=c("blue","red"),lwd=1,pch=1)' % int(int(all_junc[-1])/1000), file=OUT)
                print("dev.off()", file=OUT)

        
        def annotate_junction(self,refgene,outfile=None,min_intron=50):
                '''Annotate splicing junctions in SAM file. Note that a (long) read might have multiple splicing
                events  (splice multiple times), and the same splicing events can be consolidated into a single
                junction'''
                
                if outfile is None:
                        out_file = self.fileName + ".junction.xls"
                        out_file2 = self.fileName + ".junction_plot.r"
                else:
                        out_file = outfile + ".junction.xls"
                        out_file2 = outfile + ".junction_plot.r"
                if refgene is None:
                        print("You must provide reference gene model in bed format.", file=sys.stderr)
                        sys.exit(1)
                OUT = open(out_file,'w')
                ROUT = open(out_file2,'w')
                
                #reading reference gene model
                refIntronStarts=collections.defaultdict(dict)
                refIntronEnds=collections.defaultdict(dict)     
                total_junc =0
                novel35_junc =0
                novel3or5_junc =0
                known_junc =0
                splicing_events=collections.defaultdict(int)    
                
                print("\treading reference bed file: ",refgene, " ... ", end=' ', file=sys.stderr)
                for line in open(refgene,'r'):
                        if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                        fields = line.split()
                        if(len(fields)<12):
                                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                                continue
                        chrom     = fields[0].upper()
                        tx_start = int( fields[1] )
                        tx_end   = int( fields[2] )
                        if int(fields[9] ==1):
                                continue        
                        
                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        intron_start = exon_ends[:-1]
                        intron_end=exon_starts[1:]
                        for i_st,i_end in zip (intron_start, intron_end):
                                refIntronStarts[chrom][i_st] =i_st
                                refIntronEnds[chrom][i_end] =i_end                      
                print("Done", file=sys.stderr)
                
                #reading input SAM file
                print("\tProcessing "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        
                        if (flagCode & 0x0004) != 0: continue                           #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):            #skip multiple mapped reads
                                continue
                        if (len(ParseSAM._splicedHit_pat.findall(fields[5]))==1):       #skip non-spilced reads
                                continue

                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        blockStart=[]
                        blockSize=[]
                        blockEnd=[]             
                        #if intron size < min_intron, skip. 
                        flag=0
                        for i in range(1,len(comb),2):
                                if comb[i] < min_intron:
                                        flag=1
                                        break
                        if flag ==1:
                                continue                
                        
                        total_junc += (len(comb) -1)/2
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        for st,size in zip(blockStart,blockSize):
                                end = st + size
                                blockEnd.append(end)
                        intron_st = blockEnd[:-1]
                        intron_end = blockStart[1:]
                        for i_st,i_end in zip(intron_st, intron_end):
                                splicing_events[chrom + ":" + str(i_st) + ":" + str(i_end)] += 1
                                if (i_st in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                                        known_junc +=1                                                                                                                                          #known both
                                elif (i_st not in refIntronStarts[chrom] and i_end not in refIntronEnds[chrom]):
                                        novel35_junc +=1                                                                                                                                
                                else:
                                        novel3or5_junc +=1
                #self.f.seek(0)
                print("Done", file=sys.stderr)
                
                print('pdf("splicing_events_pie.pdf")', file=ROUT)
                print("events=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc)])+ ')', file=ROUT)
                print('pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing events",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc)), file=ROUT)
                print("dev.off()", file=ROUT)
                
                print("\n===================================================================", file=sys.stderr)
                print("Total splicing  Events:\t" + str(total_junc), file=sys.stderr)
                print("Known Splicing Events:\t" + str(known_junc), file=sys.stderr)
                print("Partial Novel Splicing Events:\t" + str(novel3or5_junc), file=sys.stderr)
                print("Novel Splicing Events:\t" + str(novel35_junc), file=sys.stderr)
                
                #reset variables
                total_junc =0
                novel35_junc =0
                novel3or5_junc =0
                known_junc =0
                
                print("chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation", file=OUT)
                for i in splicing_events:
                        total_junc += 1
                        (chrom, i_st, i_end) = i.split(":")
                        print('\t'.join([chrom.replace("CHR","chr"),i_st,i_end]) + '\t' + str(splicing_events[i]) + '\t', end=' ', file=OUT)
                        i_st = int(i_st)
                        i_end = int(i_end)
                        if (i_st in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                                print("annotated", file=OUT)
                                known_junc +=1
                        elif (i_st not in refIntronStarts[chrom] and i_end not in refIntronEnds[chrom]):
                                print('complete_novel', file=OUT)
                                novel35_junc +=1
                        else:
                                print('partial_novel', file=OUT)
                                novel3or5_junc +=1
                
                
                print("\nTotal splicing  Junctions:\t" + str(total_junc), file=sys.stderr)
                print("Known Splicing Junctions:\t" + str(known_junc), file=sys.stderr)
                print("Partial Novel Splicing Junctions:\t" + str(novel3or5_junc), file=sys.stderr)
                print("Novel Splicing Junctions:\t" + str(novel35_junc), file=sys.stderr)
                print("\n===================================================================", file=sys.stderr)
                
                print('pdf("splicing_junction_pie.pdf")', file=ROUT)
                print("junction=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc,)])+ ')', file=ROUT)
                print('pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing junctions",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc)), file=ROUT)
                print("dev.off()", file=ROUT)
                #print >>ROUT, "mat=matrix(c(events,junction),byrow=T,ncol=3)"
                #print >>ROUT, 'barplot(mat,beside=T,ylim=c(0,100),names=c("known","partial\nnovel","complete\nnovel"),legend.text=c("splicing events","splicing junction"),ylab="Percent")'

        def mRNA_RPKM(self,refbed,outfile=None):
                '''calculate mRNA's RPKM value'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if outfile is None:
                        rpkm_file = self.fileName + ".RPKM.xls"
                else:
                        rpkm_file = outfile + ".RPKM.xls"               
                RPKM_OUT = open(rpkm_file,'w')                  
                
                ranges={}
                totalReads=0
                cUR_num = 0     #number

                RPKM_table={}
                rawCount_table={}
                mRNAlen_table={}
                RPKM_head=['chr','start','end','name','score','strand','length','rawCount','RPKM']              
        
                #read SAM 
                print("Reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue           #skip unmap reads
                        totalReads +=1
                        if not ParseSAM._uniqueHit_pat.search(line):            #skip multiple mapped reads
                                continue

                        chrom = fields[2].upper()
                        chromStart = string.atoi(fields[3])-1
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        cUR_num += (len(comb) +1)/2
                        blockStart=[]
                        blockSize=[]
                        
                        for i in range(0,len(comb),2):
                                blockStart.append(chromStart + sum(comb[:i]) )
                                
                        for i in range(0,len(comb),2):
                                blockSize.append(comb[i])
                        
                        for st,size in zip(blockStart,blockSize):
                                mid = int(st) + (size/2)
                                if chrom not in ranges:
                                        ranges[chrom] = Intersecter()
                                ranges[chrom].add_interval( Interval( mid, mid ) )
                                        
                print("Done", file=sys.stderr)

                
                print("Calculating RPKM ...", end=' ', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith(('#','track','browser')):continue  
                                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                                continue
                        mRNA_count=0    #we need to initializ it to 0 for each gene
                        mRNA_len=sum(exon_sizes)
                        for st,end in zip(exon_starts,exon_ends):
                                if chrom in ranges:
                                        mRNA_count += len(ranges[chrom].find(st,end))                                   
                        mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * cUR_num)
                        
                        mRNAlen_table[key] = mRNA_len
                        RPKM_table[key] = str(mRNA_RPKM)
                        rawCount_table[key] = str(mRNA_count)
                print("Done", file=sys.stderr)
                
                print('\t'.join(RPKM_head), file=RPKM_OUT)
                for k in RPKM_table:
                        print(k + '\t', end=' ', file=RPKM_OUT)
                        print(str(mRNAlen_table[k]) + '\t', end=' ', file=RPKM_OUT)
                        print(str(rawCount_table[k]) + '\t', end=' ', file=RPKM_OUT)
                        print(str(RPKM_table[k]) + '\t', file=RPKM_OUT)
                return RPKM_table
                self.f.seek(0)  

        def strand_specificity(self,refbed,genome,outfile=None,sp="GTAG,GCAG,ATAC"):
                '''Measure the strand specificity of strand specific RNA-seq protocol. For non-splice read,
                use the parental gene as standard, for spliced read, use the splicing motif as strandard'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                if genome is None:
                        print("You must specify genome sequence in fasta format\n", file=sys.stderr)
                        exit(0)
                
                if outfile is None:
                        strand_file = self.fileName + ".strand.infor"
                else:
                        strand_file = outfile + ".strand.infor"         
                OUT = open(strand_file,'w')
                print("read_type\tread_id\tread_seq\tchr\tStart\tCigar\tprotocol_strand\tgene_strand", file=OUT)        
                
                transtab = string.maketrans("ACGTNX","TGCANX")
                motif=sp.upper().split(',')
                motif_rev = [m.translate(transtab)[::-1] for m in motif]
                
                #load genome
                print("\tloading "+genome+'...', file=sys.stderr)
                tmp=fasta.Fasta(genome)
                
                #load reference gene model
                gene_ranges={}
                print("reading reference gene model ...", end=' ', file=sys.stderr)
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
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                                continue
                        if chrom not in gene_ranges:
                                gene_ranges[chrom]=Intersecter()
                        gene_ranges[chrom].insert(tx_start,tx_end,strand)                                                       
                print("Done", file=sys.stderr)          

                #read SAM 
                
                read_type="unknown"
                strand_from_protocol = 'unknown'
                strand_from_gene='unknown'
                strand_stat=collections.defaultdict(int)
                print("Reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue                                           #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):continue            #skip multiple mapped reads
                        if (flagCode & 0x0100 !=0): continue                                            #skip non primary hit
                        if (flagCode & 0x0200 !=0): continue                                            #skip QC-failed
                        if (flagCode & 0x0400 !=0): continue                                            #skip PCR artifact
                        if (flagCode & 0x0010 !=0): strand_from_protocol = '-'
                        if (flagCode & 0x0010 ==0): strand_from_protocol = '+'
                        if (flagCode & 0x0040 !=0): read_type="read_1"
                        if (flagCode & 0x0080 !=0): read_type="read_2"
                        chrom = fields[2]
                        comb=[int(i) for i in ParseSAM._splicedHit_pat.findall(fields[5])]      #"9M4721N63M3157N8M" return ['9', '4721', '63', '3157', '8']
                        readStart = string.atoi(fields[3])-1
                        
                        
                        #for non spliced read
                        if len(comb)==1:
                                readEnd = readStart + len(fields[9])
                                if chrom in gene_ranges:
                                        if len(set(gene_ranges[chrom].find(readStart,readEnd)))>1:    
                                                strand_from_gene="overlap"
                                        elif len(set(gene_ranges[chrom].find(readStart,readEnd)))==1:
                                                strand_from_gene = set(gene_ranges[chrom].find(readStart,readEnd)).pop()
                                        else:
                                                strand_from_gene="intergenic"
                                
                                print(read_type + '\t' + fields[0] + '\t' + fields[9] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[5] +'\t', end=' ', file=OUT)
                                print(strand_from_protocol + '\t' + strand_from_gene, file=OUT)           
                                strand_stat[read_type + '\t' + strand_from_protocol +'\t' + strand_from_gene] +=1  


                        #for spliced read
                        if len(comb)>=3:
                                splice_strand=[]
                                blockStart=[]
                                blockSize=[]
                                blockEnd =[]
                                for i in range(0,len(comb),2):
                                        blockStart.append(readStart + sum(comb[:i]) )
                                for i in range(0,len(comb),2):
                                        blockSize.append(comb[i])
                                blockEnd=list(map((lambda x,y:x+y),blockStart,blockSize))
                                intron_start=blockEnd[:-1]
                                intron_end=blockStart[1:]
                                for st,end in zip(intron_start,intron_end):
                                        try:
                                                splice_motif = str(tmp.fetchSeq(chrom, st, st+2)) + str(tmp.fetchSeq(chrom, end-2,end))
                                        except:
                                                        print(line)
                                        if splice_motif in motif:
                                                splice_strand.append('+')
                                        elif splice_motif in motif_rev:
                                                splice_strand.append('-')
                                        else:
                                                splice_strand.append('unknown motif')


                                if len(set(splice_strand))>1:
                                                strand_from_splice = 'unknown motif'
                                else:
                                                strand_from_splice = set(splice_strand).pop()
                                print(read_type + '\t' + fields[0] + '\t' + fields[9] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[5] +'\t', end=' ', file=OUT)
                                print(strand_from_protocol + '\t' + strand_from_splice, file=OUT)
                                
                                strand_stat[read_type + '\t' + strand_from_protocol +'\t' + strand_from_splice] +=1
                                                        
                print("Done", file=sys.stderr)
                
                print("read_type\tstrand_expected\tstrand_observed\tcount")
                for i in sorted(strand_stat):
                                print(str(i) +'\t' + str(strand_stat[i]))
                
        def clipping_profile(self,outfile=None):
                '''calculate profile of soft clipping'''
                if outfile is None:
                        out_file1 = self.fileName + ".clipping_profile.xls"
                        out_file2 = self.fileName + ".clipping_profile.r"
                else:
                        out_file1 = outfile + ".clipping_profile.xls"
                        out_file2 = outfile + ".clipping_profile.r"
                
                OUT=open(out_file1,'w')
                ROUT=open(out_file2,'w')
                print("Position\tRead_Total\tRead_clipped", file=OUT)
                soft_p = re.compile(r'(.*?)(\d+)S')
                read_part = re.compile(r'(\d+)[MIS=X]')
                total_read =0
                skip_part_of_read =0
                soft_clip_profile=collections.defaultdict(int)
                
                read_pos=[]
                clip_count=[]
                print("Reading "+ self.fileName + '...', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue                                           #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):continue            #skip multiple mapped reads
                        if (flagCode & 0x0100 !=0): continue                                            #skip non primary hit
                        if (flagCode & 0x0200 !=0): continue                                            #skip QC-failed
                        if (flagCode & 0x0400 !=0): continue                                            #skip PCR artifact
                        total_read +=1
                        m = soft_p.findall(fields[5])
                        
                        skip_part_of_read =0
                        if len(m)==0:                                   #NO soft clip
                                continue
                        else:
                                for j in m: 
                                        skip_part_of_read += sum([int(i) for i in read_part.findall(j[0])])
                                        for n in range(skip_part_of_read,(skip_part_of_read + int(j[1]))):
                                                soft_clip_profile[n]+=1
                                        skip_part_of_read += int(j[1])
                for i in soft_clip_profile:
                        read_pos.append(str(i))
                        clip_count.append(str(soft_clip_profile[i]))
                        print(str(i) + '\t' + str(total_read) + '\t' + str(soft_clip_profile[i]), file=OUT)
                print("pdf(\"%s\")" % outfile + '.clipping_profile.pdf', file=ROUT)
                print("read_pos=c(" + ','.join(read_pos) + ')', file=ROUT)
                print("count=c(" + ','.join(clip_count) + ')', file=ROUT)
                print('plot(read_pos,1-(count/%d),col="blue",main="clipping profile",xlab="Position of reads",ylab="Mappability",type="b")' % total_read, file=ROUT)
                print("dev.off()", file=ROUT)
                
        def insertion_profile(self,read_len,outfile=None):
                '''calculate profile of insertion (insertion means insertion to the reference)'''
                if outfile is None:
                        out_file1 = self.fileName + ".insertion_profile.xls"
                        out_file2 = self.fileName + ".insertion_profile.r"
                else:
                        out_file1 = outfile + ".insertion_profile.xls"
                        out_file2 = outfile + ".insertion_profile.r"
                
                OUT=open(out_file1,'w')
                ROUT=open(out_file2,'w')
                print("Position\tRead_Total\tRead_clipped", file=OUT)
                soft_p = re.compile(r'(.*?)(\d+)I')
                read_part = re.compile(r'(\d+)[MIS=X]')
                total_read =0
                skip_part_of_read =0
                soft_clip_profile=collections.defaultdict(int)
                print("Reading "+ self.fileName + '...', end=' ', file=sys.stderr)
                for line in self.f:
                        if line.startswith("@"):continue
                        fields=line.rstrip('\n ').split()
                        flagCode=string.atoi(fields[1])
                        if (flagCode & 0x0004) != 0: continue                                           #skip unmap reads
                        if not ParseSAM._uniqueHit_pat.search(line):continue            #skip multiple mapped reads
                        if (flagCode & 0x0100 !=0): continue                                            #skip non primary hit
                        if (flagCode & 0x0200 !=0): continue                                            #skip QC-failed
                        if (flagCode & 0x0400 !=0): continue                                            #skip PCR artifact
                        total_read +=1
                        m = soft_p.findall(fields[5])
                        
                        skip_part_of_read =0
                        if len(m)==0:                                   #NO soft clip
                                continue
                        else:
                                for j in m: 
                                        skip_part_of_read += sum([int(i) for i in read_part.findall(j[0])])
                                        for n in range(skip_part_of_read,(skip_part_of_read + int(j[1]))):
                                                soft_clip_profile[n]+=1
                                        skip_part_of_read += int(j[1])
                for i in range(0,read_len):
                        print(str(i) + '\t' + str(total_read) + '\t' + str(soft_clip_profile[i]), file=OUT)

class ParseBAM(object):
        '''This class provides fuctions to parsing/processing/transforming SAM or BAM files. The input
        file could be either SAM or BAM format file'''
        
        multi_hit_tags=['H0','H1','H2','IH','NH']
        def __init__(self,inputFile):
                '''constructor. input could be bam or sam'''
                try:
                        self.samfile = pysam.Samfile(inputFile,'rb')
                        if len(self.samfile.header) ==0:
                                print("BAM/SAM file has no header section. Exit!", file=sys.stderr)
                                sys.exit(1)
                        self.bam_format = True
                except:
                        self.samfile = pysam.Samfile(inputFile,'r')
                        if len(self.samfile.header) ==0:
                                print("BAM/SAM file has no header section. Exit!", file=sys.stderr)
                                sys.exit(1)
                        self.bam_format = False

        def stat (self, q_cut=30):
                '''Calculate mapping statistics'''
                R_total=0
                R_qc_fail=0
                R_duplicate=0
                R_nonprimary=0
                R_unmap =0
                
                R_multipleHit=0
                R_uniqHit=0     #all the following count should be sum to uniqHit
                
                R_read1=0
                R_read2=0
                R_reverse =0
                R_forward=0
                R_nonSplice=0
                R_splice=0
                R_properPair =0 
                R_pair_diff_chrom = 0
                R_mitochondria = 0

                if self.bam_format:
                        print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:
                        print("Load SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                flag=0
                                aligned_read = next(self.samfile)
                                R_total +=1
                                if aligned_read.is_qcfail:                      #skip QC fail read
                                        R_qc_fail +=1
                                        continue
                                if aligned_read.is_duplicate:           #skip duplicate read
                                        R_duplicate +=1
                                        continue
                                if aligned_read.is_secondary:           #skip non primary hit
                                        R_nonprimary +=1
                                        continue
                                if aligned_read.is_unmapped:            #skip unmap read
                                        R_unmap +=1
                                        continue                
                                        
                                if aligned_read.mapq < q_cut:
                                        R_multipleHit +=1
                                        continue                                                #skip multiple map read                         
                                if aligned_read.mapq >= q_cut:
                                        R_uniqHit +=1
                                        
                                        if aligned_read.is_read1:
                                                R_read1 +=1
                                        if aligned_read.is_read2:
                                                R_read2 +=1
                                        if aligned_read.is_reverse:
                                                R_reverse +=1
                                        else:
                                                R_forward +=1
                                        introns = bam_cigar.fetch_intron('chr1', aligned_read.pos, aligned_read.cigar)
                                        if len(introns)==0:
                                                R_nonSplice +=1
                                        if len(introns)>=1:
                                                R_splice +=1
                                        if aligned_read.is_proper_pair:
                                                R_properPair +=1
                                                R_read1_ref = self.samfile.getrname(aligned_read.tid)
                                                R_read2_ref = self.samfile.getrname(aligned_read.rnext)
                                                if R_read1_ref != R_read2_ref:
                                                        #print aligned_read.qname
                                                        R_pair_diff_chrom +=1
                                
                except StopIteration:
                        print("Done", file=sys.stderr)          
                #self.samfile.seek(current_pos)
                                
                print("\n#==================================================", file=sys.stdout)
                print("#All numbers are READ count", file=sys.stdout)
                print("#==================================================\n", file=sys.stdout)
                print("%-40s%d" % ("Total records:",R_total), file=sys.stdout)
                print("\n", end=' ', file=sys.stdout)
                print("%-40s%d" % ("QC failed:",R_qc_fail), file=sys.stdout)
                print("%-40s%d" % ("Optical/PCR duplicate:", R_duplicate), file=sys.stdout)
                print("%-40s%d" % ("Non primary hits", R_nonprimary), file=sys.stdout)
                print("%-40s%d" % ("Unmapped reads:",R_unmap), file=sys.stdout)
                print("%-40s%d" % ("mapq < mapq_cut (non-unique):",R_multipleHit), file=sys.stdout)
                print("\n", end=' ', file=sys.stdout)
                print("%-40s%d" % ("mapq >= mapq_cut (unique):",R_uniqHit), file=sys.stdout)
                print("%-40s%d" % ("Read-1:",R_read1), file=sys.stdout)
                print("%-40s%d" % ("Read-2:",R_read2), file=sys.stdout)
                print("%-40s%d" % ("Reads map to '+':",R_forward), file=sys.stdout)
                print("%-40s%d" % ("Reads map to '-':",R_reverse), file=sys.stdout)
                print("%-40s%d" % ("Non-splice reads:",R_nonSplice), file=sys.stdout)
                print("%-40s%d" % ("Splice reads:",R_splice), file=sys.stdout)  
                print("%-40s%d" % ("Reads mapped in proper pairs:",R_properPair), file=sys.stdout)
                print("%-40s%d" % ("Proper-paired reads map to different chrom:",R_pair_diff_chrom), file=sys.stdout)
        
        def configure_experiment(self,refbed,sample_size, q_cut = 30):
                '''Given a BAM/SAM file, this function will try to guess the RNA-seq experiment:
                        1) single-end or pair-end
                        2) strand_specific or not
                        3) if it is strand-specific, what's the strand_ness of the protocol
                '''
                
                        #how many reads you want to sample
                count =0
                p_strandness=collections.defaultdict(int)
                s_strandness=collections.defaultdict(int)
                #load reference gene model
                gene_ranges={}
                print("Reading reference gene model " + refbed + ' ...', end=' ', file=sys.stderr)
                for line in open(refbed,'r'):
                        try:
                                if line.startswith(('#','track','browser')):continue  
                                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0]
                                txStart  = int( fields[1] )
                                txEnd    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5]
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                                continue
                        if chrom not in gene_ranges:
                                gene_ranges[chrom]=Intersecter()
                        gene_ranges[chrom].insert(txStart,txEnd,strand)                                                 
                print("Done", file=sys.stderr)          
                
                #read SAM/BAM file
                #current_pos = self.samfile.tell()
                print("Loading SAM/BAM file ... ", end=' ', file=sys.stderr)
                try:
                        while(1):
                                if count >= sample_size:
                                        break
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:                      #skip low quanlity
                                        continue
                                if aligned_read.is_duplicate:           #skip duplicate read
                                        continue
                                if aligned_read.is_secondary:           #skip non primary hit
                                        continue
                                if aligned_read.is_unmapped:            #skip unmap read
                                        continue                
                                if aligned_read.mapq < q_cut:
                                        continue                                                                                                                
                                
                                chrom = self.samfile.getrname(aligned_read.tid)
                                if aligned_read.is_paired:
                                        if aligned_read.is_read1:
                                                read_id = '1'
                                        if aligned_read.is_read2:
                                                read_id = '2'
                                        if aligned_read.is_reverse:
                                                map_strand = '-'
                                        else:
                                                map_strand = '+'
                                        readStart = aligned_read.pos
                                        readEnd = readStart + aligned_read.qlen
                                        if chrom in gene_ranges:
                                                tmp = set(gene_ranges[chrom].find(readStart,readEnd))
                                                if len(tmp) == 0: continue
                                                strand_from_gene = ':'.join(tmp)
                                                p_strandness[read_id + map_strand + strand_from_gene]+=1        
                                                count += 1
                                else:
                                        if aligned_read.is_reverse:
                                                map_strand = '-'
                                        else:
                                                map_strand = '+'                                        
                                        readStart = aligned_read.pos
                                        readEnd = readStart + aligned_read.qlen
                                        if chrom in gene_ranges:
                                                tmp = set(gene_ranges[chrom].find(readStart,readEnd))
                                                if len(tmp) == 0: continue
                                                strand_from_gene = ':'.join(tmp)
                                                s_strandness[map_strand + strand_from_gene]+=1
                                                count += 1

                except StopIteration:
                        print("Finished", file=sys.stderr)              
                #self.samfile.seek(current_pos)
                
                print("Total " + str(count) + " usable reads were sampled", file=sys.stderr)
                protocol="unknown"
                strandness=None
                spec1=0.0
                spec2=0.0
                other=0.0
                if len(p_strandness) >0 and len(s_strandness) ==0 :
                        protocol="PairEnd"
                        #for k,v in p_strandness.items():
                        #       print >>sys.stderr, k + '\t' + str(v)
                        spec1= (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
                        spec2= (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
                        other = 1-spec1-spec2
                        
                elif len(s_strandness) >0 and len(p_strandness) ==0 :
                        protocol="SingleEnd"
                        #for k,v in s_strandness.items():
                        #       print  >>sys.stderr, k + '\t' + str(v)
                        spec1 = (s_strandness['++'] + s_strandness['--'])/float(sum(s_strandness.values()))
                        spec2 = (s_strandness['+-'] + s_strandness['-+'])/float(sum(s_strandness.values()))
                        other = 1-spec1-spec2
                else:
                        protocol="Mixture"
                        spec1 = "NA"
                        spec2 = "NA"
                        other = "NA"
                return [protocol,spec1,spec2,other]

        def bamTowig(self,outfile,chrom_sizes, chrom_file,skip_multi=True,strand_rule=None,WigSumFactor=None,q_cut=30):
                """Convert BAM/SAM file to wig file. chrom_size is dict with chrom as key and chrom_size as value
                strandRule should be determined from \"infer_experiment\". such as \"1++,1--,2+-,2-+\". When
                WigSumFactor is provided, output wig file will be normalized to this number """
                
                #strand_rule={'1+':'-','1-':'+','2+':'+','2-':'-'}
                strandRule={}
                if strand_rule is None:                                                                                                 # Not strand-specific
                        pass                                                                                                                            
                elif len(strand_rule.split(',')) ==4:                                                                   #PairEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
                elif len(strand_rule.split(',')) ==2:                                                                   #singeEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]]=i[1]
                else:
                        print("Unknown value of option :'strand_rule' " + strand_rule, file=sys.stderr)
                        sys.exit(1)
                if len(strandRule) == 0:
                        FWO = open(outfile + '.wig','w')
                else:
                        FWO = open(outfile + '.Forward.wig','w')
                        RVO = open(outfile + '.Reverse.wig','w')
                
                read_id=''
                
                for chr_name, chr_size in list(chrom_sizes.items()):            #iterate each chrom
                        try:
                                self.samfile.fetch(chr_name,0,chr_size)
                        except:
                                print("No alignments for " + chr_name + '. skipped', file=sys.stderr)
                                continue
                        print("Processing " + chr_name + " ...", file=sys.stderr)
                        if len(strandRule) == 0: FWO.write('variableStep chrom='+chr_name+'\n')
                        else:
                                FWO.write('variableStep chrom='+chr_name+'\n')
                                RVO.write('variableStep chrom='+chr_name+'\n')
                        Fwig = collections.defaultdict(int)
                        Rwig = collections.defaultdict(int)
                        alignedReads = self.samfile.fetch(chr_name,0,chr_size)          
                        for aligned_read in alignedReads:
                                if aligned_read.is_qcfail:continue                      #skip low quanlity
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                
                                if skip_multi:
                                        if aligned_read.mapq < q_cut:
                                                continue                                
                                if aligned_read.is_paired:
                                        if aligned_read.is_read1:read_id = '1'
                                        if aligned_read.is_read2:read_id = '2'
                                
                                if aligned_read.is_reverse:map_strand = '-'
                                else:map_strand = '+'
                                
                                key = read_id + map_strand
                                
                                hit_st = aligned_read.pos
                                for block in bam_cigar.fetch_exon(chr_name, hit_st, aligned_read.cigar): 
                                        for pos in range(block[1]+1,block[2]+1):        
                                                if len(strandRule) == 0: Fwig[pos] +=1.0        #this is NOT strand specific. everything into Fwig
                                                else:                                                                           #this is strand specific. separate Fwig and Rwig
                                                        if strandRule[key] == '+':Fwig[pos] +=1.0
                                                        if strandRule[key] == '-':Rwig[pos] -=1.0
                        if WigSumFactor is None:        #not normalize  
                                if len(strandRule) == 0:                                                        #this is NOT strand specific.
                                        for pos in sorted (Fwig.keys()):
                                                print("%d\t%.2f" % (pos,Fwig[pos]), file=FWO)
                                else:
                                        for pos in sorted (Fwig.keys()):
                                                print("%d\t%.2f" % (pos,Fwig[pos]), file=FWO)
                                        for pos in sorted (Rwig.keys()):
                                                print("%d\t%.2f" % (pos,Rwig[pos]), file=RVO)
                        else:                                           #normalize wig signal to WigSumFactor
                                if len(strandRule) == 0:                                                        #this is NOT strand specific.
                                        for pos in sorted (Fwig.keys()):
                                                print("%d\t%.2f" % (pos,Fwig[pos]*WigSumFactor), file=FWO)
                                else:
                                        for pos in sorted (Fwig.keys()):
                                                print("%d\t%.2f" % (pos,Fwig[pos]*WigSumFactor), file=FWO)
                                        for pos in sorted (Rwig.keys()):
                                                print("%d\t%.2f" % (pos,Rwig[pos]*WigSumFactor), file=RVO)
                FWO.close()
                if len(strandRule) != 0:
                        RVO.close()
                if len(strandRule) == 0:
                        try:
                                import subprocess
                                print("Run " + "wigToBigWig " + outfile + '.wig ' + chrom_file + ' ' +  outfile + ".bw ")
                                subprocess.call("wigToBigWig -clip " + outfile + '.wig ' + chrom_file + ' ' +  outfile + ".bw ",shell=True)
                        except:
                                print("Failed to call \"wigToBigWig\".", file=sys.stderr)
                                pass
                else:
                        try:
                                import subprocess
                                subprocess.call("wigToBigWig -clip " + outfile + '.Forward.wig ' + chrom_file + ' ' +  outfile + ".Forward.bw ",shell=True)
                                subprocess.call("wigToBigWig -clip " + outfile + '.Reverse.wig ' + chrom_file + ' ' +  outfile + ".Reverse.bw ",shell=True)
                        except:
                                print("Failed to call \"wigToBigWig\".", file=sys.stderr)
                                pass                    

        def calWigSum(self,chrom_sizes, skip_multi=True):
                """Calculate wigsum from BAM file"""
                
                print("Calcualte wigsum ... ", file=sys.stderr)
                wigsum = 0.0
                read_id=''
                for chr_name, chr_size in list(chrom_sizes.items()):            #iterate each chrom
                        try:
                                self.samfile.fetch(chr_name,0,chr_size)
                        except:
                                print("No alignments for " + chr_name + '. skipped', file=sys.stderr)
                                continue
                        print("Processing " + chr_name + " ...", file=sys.stderr)

                        alignedReads = self.samfile.fetch(chr_name,0,chr_size)          
                        for aligned_read in alignedReads:
                                flag=0
                                if aligned_read.is_qcfail:continue                      #skip low quanlity
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                
                                if skip_multi:
                                        if len(aligned_read.tags)>0:            #( ("NM", 1),("RG", "L1") )
                                                for i in aligned_read.tags:
                                                        if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
                                                                flag=1                                          #multiple hit read
                                                                break
                                        if flag==1:continue                                             #skip multiple map read         
                                
                                if aligned_read.is_paired:
                                        if aligned_read.is_read1:read_id = '1'
                                        if aligned_read.is_read2:read_id = '2'
                                
                                if aligned_read.is_reverse:map_strand = '-'
                                else:map_strand = '+'
                                
                                key = read_id + map_strand
                                
                                hit_st = aligned_read.pos
                                for block in bam_cigar.fetch_exon(chr_name, hit_st, aligned_read.cigar): 
                                        wigsum += (block[2] - block[1])
                return wigsum   

        def bam2fq(self,prefix, paired = True):
                """Convert BAM/SAM into fastq files"""
                
                transtab = string.maketrans("ACGTNX","TGCANX")
                
                if paired:
                        OUT1 = open(prefix + '.R1.fastq','w')
                        OUT2 = open(prefix + '.R2.fastq','w')
                        read1_count = 0
                        read2_count = 0
                else:
                        OUT = open(prefix + '.fastq','w')
                        read_count = 0
                read_name = ''
                read_seq = ''
                read_qual = ''

                print("Convert BAM/SAM file into fastq format ... ", end=' ', file=sys.stderr)
                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                read_name = aligned_read.qname
                                read_seq = aligned_read.seq.upper()
                                read_qual = aligned_read.qual
                                if aligned_read.is_reverse:
                                        read_seq = read_seq.translate(transtab)[::-1]
                                        read_qual = read_qual[::-1]
                                if paired:
                                        if aligned_read.is_read1:
                                                read1_count += 1
                                                if not read_name.endswith('/1'): 
                                                        print('@' + read_name + '/1', file=OUT1)
                                                else:
                                                        print('@' + read_name, file=OUT1)
                                                print(read_seq, file=OUT1)
                                                print('+', file=OUT1)
                                                print(read_qual, file=OUT1)
                                        if aligned_read.is_read2:
                                                read2_count += 1
                                                if not read_name.endswith('/2'): 
                                                        print('@' + read_name + '/2', file=OUT2)
                                                else:
                                                        print('@' + read_name, file=OUT2)
                                                print(read_seq, file=OUT2)
                                                print('+', file=OUT2)
                                                print(read_qual, file=OUT2)     
                                else:           #single end
                                        read_count += 1
                                        print('@' + read_name, file=OUT)
                                        print(read_seq, file=OUT)
                                        print('+', file=OUT)
                                        print(read_qual, file=OUT)      
                                        
                except StopIteration:
                        print("Done", file=sys.stderr)
                if paired: 
                        print("read_1 count: %d" %  read1_count, file=sys.stderr)
                        print("read_2 count: %d" %  read2_count, file=sys.stderr)
                else:
                        print("read count: %d"  % read_count, file=sys.stderr)
                                        
        def calculate_rpkm(self,geneFile,outfile,strand_rule=None):
                '''calculate RPKM vaues. For single end RNA-seq, if it is strand specific, we assume that
                read plus mapped indicates a gene on plus strand.(similar to minus). 
                Advantages: works for both SAM and BAM
                                        works for both sorted and unsorted BAM/SAM file
                                        works for both index or unindexed BAM/SAM file
                                        much faster than indexing bam file
                Disadvantage: random access BAM file was disabled, thus large mount of RAM is required
                
                strand_rule: could be the following values:
                        '1++,1--,2+-,2-+
                        '1+-,1-+,2++,2--
                        '++,--'
                        '+-,-+'
                        None
                '''
                
                strandRule={}
                if strand_rule is None:                                                                                                 # Not strand-specific
                        pass                                                                                                                            
                elif len(strand_rule.split(',')) ==4:                                                                   #PairEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
                elif len(strand_rule.split(',')) ==2:                                                                   #singeEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]]=i[1]
                else:
                        print("Unknown value of option :'strand_rule' " + strand_rule, file=sys.stderr)
                        sys.exit(1)
                
                uniq_read=0
                total_tags=0
                plus_ranges={}
                minus_ranges={}
                unstrand_ranges={}
                
                rpkm_value={}
                
                RPKM_OUT = open(outfile,'w')
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)
                
                #current_pos = self.samfile.tell()
                try:
                        while(1):
                                flag=0
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read

                                if len(aligned_read.tags)>0:            #( ("NM", 1),("RG", "L1") )
                                        for i in aligned_read.tags:
                                                if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
                                                        flag=1                                          #multiple hit read
                                                        break
                                if flag==1:continue                                             #skip multiple map read         
                                
                                uniq_read +=1
                                
                                if aligned_read.is_paired:
                                        if aligned_read.is_read1:read_id = '1'
                                        if aligned_read.is_read2:read_id = '2'
                                else:
                                        read_id = ''
                                if aligned_read.is_reverse:map_strand = '-'
                                else:map_strand = '+'
                                
                                strand_key = read_id + map_strand
                                
                                chrom = self.samfile.getrname(aligned_read.tid).upper()
                                hit_st = aligned_read.pos
                                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)
                                total_tags += len(exon_blocks)
                                                
                                #construct bitset
                                if strand_rule is not None:     
                                        if strandRule[strand_key] == '+':
                                                for block in exon_blocks:
                                                        mid = block[1] + int((block[2] - block[1])/2)
                                                        if chrom not in plus_ranges:plus_ranges[chrom] = Intersecter()
                                                        plus_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
                                        elif strandRule[strand_key] == '-':
                                                for block in exon_blocks:
                                                        mid = block[1] + int((block[2] - block[1])/2)
                                                        if chrom not in minus_ranges:minus_ranges[chrom] = Intersecter()        
                                                        minus_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
                                elif strand_rule is None:       
                                        for block in exon_blocks:
                                                mid = block[1] + int((block[2] - block[1])/2)
                                                if chrom not in unstrand_ranges:unstrand_ranges[chrom] = Intersecter()
                                                unstrand_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
                                        
                except StopIteration:
                        print("Done", file=sys.stderr)
                #self.samfile.seek(current_pos)
                print("#Total uniquely mapped reads = " + str(uniq_read), file=RPKM_OUT)
                print("#Total fragments = " + str(total_tags), file=RPKM_OUT)
                print("Assign reads to "+ geneFile + '...', end=' ', file=sys.stderr)
                for line in open(geneFile,'r'):
                        try:
                                if line.startswith('#'):continue
                                if line.startswith('track'):continue
                                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5].replace(" ","_")
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                intron_starts = exon_ends[:-1]
                                intron_ends=exon_starts[1:]
                                key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue

                                        
                        mRNA_count=0
                        mRNA_len=sum(exon_sizes)
                                
                        if (strand_rule is not None) and (strand == '-'):
                                intronNum=len(intron_starts)
                                exonNum=len(exon_starts)
                                
                                # assign reads to intron        
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in minus_ranges:
                                                hits= len(minus_ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
                                                intronNum -= 1
                                # assign reads to exon                          
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in minus_ranges:
                                                hits= len(minus_ranges[chrom].find(st,end))                                     
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
                                                exonNum -= 1
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                        elif (strand_rule is not None) and (strand == '+'):
                                intronNum=1
                                exonNum=1
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in plus_ranges:
                                                hits= len(plus_ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
                                                intronNum += 1    
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in plus_ranges:
                                                hits= len(plus_ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
                                                exonNum += 1            
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                        elif strand_rule is None:
                                intronNum=1
                                exonNum=1
                                for st,end in zip(intron_starts,intron_ends):
                                        if chrom in unstrand_ranges:
                                                hits= len(unstrand_ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
                                                intronNum += 1    
                                for st,end in zip(exon_starts,exon_ends):
                                        if chrom in unstrand_ranges:
                                                hits= len(unstrand_ranges[chrom].find(st,end))
                                                RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
                                                exonNum += 1            
                                                mRNA_count += hits
                                try:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
                                except:
                                        RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
                print("Done", file=sys.stderr)

        def readsNVC(self,outfile=None,nx=True, q_cut = 30):
                '''for each read, calculate nucleotide frequency vs position'''
                if outfile is None:
                        outfile1 = self.fileName + ".NVC.xls"
                        outfile2 = self.fileName +".NVC_plot.r"
                else:
                        outfile1 = outfile + ".NVC.xls"
                        outfile2 = outfile +".NVC_plot.r"
                FO=open(outfile1,'w')
                RS=open(outfile2,'w')
                PPcount=0
                
                transtab = string.maketrans("ACGTNX","TGCANX")
                base_freq=collections.defaultdict(int)
                a_count=[]
                c_count=[]
                g_count=[]
                t_count=[]
                n_count=[]
                x_count=[]
                if self.bam_format:print("Read BAM file ... ", end=' ', file=sys.stderr)
                else:print("Read SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                if aligned_read.mapq < q_cut: continue
                                #if aligned_read.is_unmapped:continue   #skip unmapped read
                                #if aligned_read.is_qcfail:continue     #skip low quality

                                RNA_read = aligned_read.seq.upper()             
                                if aligned_read.is_reverse:
                                        RNA_read = RNA_read.translate(transtab)[::-1]
                                for i,j in enumerate(RNA_read):
                                        key = str(i) + j
                                        base_freq[key] += 1
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                print("generating data matrix ...", file=sys.stderr)
                print("Position\tA\tC\tG\tT\tN\tX", file=FO)
                for i in range(len(RNA_read)):
                        print(str(i) + '\t', end=' ', file=FO)
                        print(str(base_freq[str(i) + "A"]) + '\t', end=' ', file=FO)
                        a_count.append(str(base_freq[str(i) + "A"]))
                        print(str(base_freq[str(i) + "C"]) + '\t', end=' ', file=FO)
                        c_count.append(str(base_freq[str(i) + "C"]))
                        print(str(base_freq[str(i) + "G"]) + '\t', end=' ', file=FO)
                        g_count.append(str(base_freq[str(i) + "G"]))
                        print(str(base_freq[str(i) + "T"]) + '\t', end=' ', file=FO)
                        t_count.append(str(base_freq[str(i) + "T"]))
                        print(str(base_freq[str(i) + "N"]) + '\t', end=' ', file=FO)
                        n_count.append(str(base_freq[str(i) + "N"]))
                        print(str(base_freq[str(i) + "X"]) + '\t', file=FO)
                        x_count.append(str(base_freq[str(i) + "X"]))
                FO.close()
                
                #generating R scripts
                print("generating R script  ...", file=sys.stderr)
                print("position=c(" + ','.join([str(i) for i in range(len(RNA_read))]) + ')', file=RS)
                print("A_count=c(" + ','.join(a_count) + ')', file=RS)
                print("C_count=c(" + ','.join(c_count) + ')', file=RS)
                print("G_count=c(" + ','.join(g_count) + ')', file=RS)
                print("T_count=c(" + ','.join(t_count) + ')', file=RS)
                print("N_count=c(" + ','.join(n_count) + ')', file=RS)
                print("X_count=c(" + ','.join(x_count) + ')', file=RS)
                
                if nx:
                        print("total= A_count + C_count + G_count + T_count + N_count + X_count", file=RS)
                        print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total) + 0.05", file=RS)
                        print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total)", file=RS)
                        
                        print('pdf(\"%s\")' % (outfile +".NVC_plot.pdf"), file=RS)
                        print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
                        print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
                        print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
                        print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
                        print('lines(position,N_count/total,type="o",pch=20,col="black")', file=RS)             
                        print('lines(position,X_count/total,type="o",pch=20,col="grey")', file=RS)      
                        print('legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C","N","X"),col=c("dark green","red","blue","cyan","black","grey"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan","black","grey"))', file=RS)
                        print("dev.off()", file=RS)
                else:
                        print("total= A_count + C_count + G_count + T_count", file=RS)
                        print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05", file=RS)
                        print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total)", file=RS)
                
                        print('pdf(\"%s\")' % (outfile +".NVC_plot.pdf"), file=RS)
                        print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
                        print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
                        print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
                        print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
                        print('legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))', file=RS)
                        print("dev.off()", file=RS)
                
                RS.close()
                #self.f.seek(0)
        def readsQual_boxplot(self,outfile,shrink=1000, q_cut=30):
                '''calculate phred quality score for each base in read (5->3)'''

                output = outfile + ".qual.r"
                FO=open(output,'w')

                if self.bam_format:print("Read BAM file ... ", end=' ', file=sys.stderr)
                else:print("Read SAM file ... ", end=' ', file=sys.stderr)

                quality = collections.defaultdict(dict) #read_pos=>quality score=>count
                q_max = -1
                q_min = 10000
                q_list=[]
                i_box={}        #key is read postion,value is 
                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                if aligned_read.mapq < q_cut: continue
                                #if aligned_read.is_unmapped:continue   #skip unmapped read
                                #if aligned_read.is_qcfail:continue             #skip low quality
                                
                                qual_str = aligned_read.qqual
                                read_len = aligned_read.rlen
                                if aligned_read.is_reverse:
                                        qual_str = qual_str[::-1]

                                for i,j in enumerate(qual_str):
                                        q=ord(j)-33
                                        if q > q_max: q_max = q
                                        if q < q_min: q_min = q
                                        try:
                                                quality[i][q] += 1
                                        except:
                                                quality[i][q] = 1
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                for p in range(0,read_len):
                        #print str(p) + ':',
                        val=[]
                        occurrence=[]
                        for q in range(q_min,q_max+1):
                                if p in quality and q in quality[p]:
                                        val.append(str(q))                              
                                        occurrence.append(str(quality[p][q]))   
                                        q_list.append(str(quality[p][q]))
                                else:
                                        q_list.append(str(0))
                        i_box[p] = 'rep(c(' + ','.join(val) + '),times=c(' + ','.join(occurrence) + ')/' + str(shrink)+ ')'
                
                
                #generate R script for boxplot
                print("pdf(\'%s\')" % (outfile + ".qual.boxplot.pdf"), file=FO)
                for i in sorted(i_box):
                        print('p'+str(i) + '<-' + i_box[i], file=FO)
                print('boxplot(' + ','.join(['p'+str(i) for i in i_box]) + ',xlab=\"Position of Read(5\'->3\')\",ylab=\"Phred Quality Score\",outline=F' + ')', file=FO)
                print("dev.off()", file=FO)
                
                
                #generate R script for heatmap
                print('\n', file=FO)
                print("pdf(\'%s\')" % (outfile + ".qual.heatmap.pdf"), file=FO)
                print("qual=c(" + ','.join(q_list)  + ')', file=FO)
                print("mat=matrix(qual,ncol=%s,byrow=F)" % (read_len), file=FO)
                print('Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c(\'spline\'))', file=FO)
                print("heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from=%s,to=%s),col = Lab.palette(256),scale=\"none\" )" % (q_min,q_max), file=FO)
                print('dev.off()', file=FO)
                
        def readGC(self,outfile=None, q_cut=30):
                '''GC content distribution of reads'''
                if outfile is None:
                        outfile1 = self.fileName + ".GC.xls"
                        outfile2 = self.fileName +".GC_plot.r"
                else:
                        outfile1 = outfile + ".GC.xls"
                        outfile2 = outfile + ".GC_plot.r"
                FO=open(outfile1,'w')
                RS=open(outfile2,'w')
                
                gc_hist=collections.defaultdict(int)    #key is GC percent, value is count of reads

                if self.bam_format:print("Read BAM file ... ", end=' ', file=sys.stderr)
                else:print("Read SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                if aligned_read.is_unmapped:continue    #skip unmapped read
                                if aligned_read.is_qcfail:continue              #skip low quality
                                if aligned_read.mapq < q_cut: continue
                                RNA_read = aligned_read.seq.upper()             
                                gc_percent = "%4.2f" % ((RNA_read.count('C') + RNA_read.count('G'))/(len(RNA_read)+0.0)*100)
                                #print gc_percent
                                gc_hist[gc_percent] += 1
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                print("writing GC content ...", file=sys.stderr)        
                print("GC%\tread_count", file=FO)
                for i in list(gc_hist.keys()):
                        print(i + '\t' + str(gc_hist[i]), file=FO)
                        
                print("writing R script ...", file=sys.stderr)
                print("pdf(\"%s\")" % (outfile +  ".GC_plot.pdf"), file=RS)
                print('gc=rep(c(' + ','.join([i for i in list(gc_hist.keys())]) + '),' + 'times=c(' + ','.join([str(i) for i in list(gc_hist.values())]) + '))', file=RS)
                print('hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100, file=RS)
                #print >>RS, "lines(density(gc),col='red')"
                print("dev.off()", file=RS)             
                #self.f.seek(0)
        def readDupRate(self,q_cut, outfile=None,up_bound=500):
                '''Calculate reads's duplicate rates'''
                if outfile is None:
                        outfile1 = self.fileName + ".seq.DupRate.xls"
                        outfile2 = self.fileName + ".pos.DupRate.xls"
                        outfile3 = self.fileName + ".DupRate_plot.r"
                else:
                        outfile1 = outfile + ".seq.DupRate.xls"
                        outfile2 = outfile + ".pos.DupRate.xls"
                        outfile3 = outfile +".DupRate_plot.r"
                SEQ=open(outfile1,'w')
                POS=open(outfile2,'w')
                RS=open(outfile3,'w')
                
                seqDup=collections.defaultdict(int)
                posDup=collections.defaultdict(int)
                
                seqDup_count=collections.defaultdict(int)
                posDup_count=collections.defaultdict(int)

                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                exon_boundary=""
                                aligned_read = next(self.samfile)
                                if aligned_read.is_unmapped:continue    #skip unmapped read
                                if aligned_read.is_qcfail:continue              #skip low quality
                                if aligned_read.mapq < q_cut: continue
                                RNA_read = aligned_read.seq.upper()             
                                seqDup[RNA_read] +=1                                    #key is read sequence

                                chrom = self.samfile.getrname(aligned_read.tid)
                                hit_st = aligned_read.pos
                                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)                           
                                for ex in exon_blocks:
                                        exon_boundary += str(ex[1]) + '-' + str(ex[2]) + ":"
                                key = chrom + ":" + str(hit_st) + ":" + exon_boundary
                                posDup[key] +=1

                except StopIteration:
                        print("Done", file=sys.stderr)

                print("report duplicte rate based on sequence ...", file=sys.stderr)
                print("Occurrence\tUniqReadNumber", file=SEQ)
                for i in list(seqDup.values()):                 #key is occurence, value is uniq reads number (based on seq)
                        seqDup_count[i] +=1
                for k in sorted(seqDup_count.keys()):   
                        print(str(k) +'\t'+ str(seqDup_count[k]), file=SEQ)
                SEQ.close()
                
                print("report duplicte rate based on mapping  ...", file=sys.stderr)
                print("Occurrence\tUniqReadNumber", file=POS)
                for i in list(posDup.values()):                 #key is occurence, value is uniq reads number (based on coord)
                        posDup_count[i] +=1
                for k in sorted(posDup_count.keys()):   
                        print(str(k) +'\t'+ str(posDup_count[k]), file=POS)
                POS.close()
                
                
                print("generate R script ...", file=sys.stderr)
                print("pdf(\'%s\')" % (outfile +".DupRate_plot.pdf"), file=RS)
                print("par(mar=c(5,4,4,5),las=0)", file=RS)
                print("seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.keys()) ]) + ')', file=RS)
                print("seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.keys()) ]) + ')', file=RS)
                print("pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.keys()) ]) + ')', file=RS)
                print("pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.keys()) ]) + ')', file=RS)
                print("plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Occurrence of read',pch=4,cex=0.8,col='blue',xlim=c(1,%d),yaxt='n')" % up_bound, file=RS)
                print("points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')", file=RS)
                print('ym=floor(max(log10(pos_uniqRead)))', file=RS)
                print("legend(%d,ym,legend=c('Sequence-based','Mapping-based'),col=c('blue','red'),pch=c(4,20))" % max(up_bound-200,1), file=RS)
                print('axis(side=2,at=0:ym,labels=0:ym)', file=RS)
                print('axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[2]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[3]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[4]*100/sum(pos_uniqRead*pos_occ))))', file=RS)
                print('mtext(4, text = "Reads %", line = 2)', file=RS)
                print('dev.off()', file=RS)
                #self.f.seek(0)
        def clipping_profile(self,outfile, q_cut, PE, type="S"):
                '''calculate profile of soft clipping or insertion'''
                
                out_file1 = outfile + ".clipping_profile.xls"
                out_file2 = outfile + ".clipping_profile.r"
                OUT = open(out_file1,'w')
                ROUT = open(out_file2,'w')
                
                print("Position\tClipped_nt\tNon_clipped_nt", file=OUT)
                
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)
                
                cigar_str = ""
                
                #single end sequencing
                if PE is False:
                        total_read = 0.0                
                        soft_clip_profile = collections.defaultdict(int)
                        try:
                                while(1):
                                        aligned_read = next(self.samfile)
                                        if aligned_read.mapq < q_cut: continue
                                        if aligned_read.is_unmapped:continue    #skip unmapped read
                                        if aligned_read.is_qcfail:continue              #skip low quality
                                
                                        total_read +=1
                                        cigar_str = bam_cigar.list2longstr(aligned_read.cigar)  # ([(0, 9), (4, 1)] ==> MMMMMMMMMS
                                        
                                        if type not in cigar_str:       # no clipping
                                                continue
                                        if aligned_read.is_reverse:
                                                cigar_str = cigar_str[::-1]             

                                        for indx,symbl in enumerate(cigar_str):
                                                if symbl == type:
                                                        soft_clip_profile[indx] += 1.0                  
                        except StopIteration:
                                print("Done", file=sys.stderr)
                        
                        print("Totoal reads used: %d" % int(total_read), file=sys.stderr)
                        read_pos = list(range(0,len(cigar_str)))
                        clip_count = []
                        for i in read_pos:
                                print(str(i) + '\t' + str(soft_clip_profile[i]) + '\t' + str(total_read - soft_clip_profile[i]), file=OUT)
                                clip_count.append(soft_clip_profile[i])
                        
                        
                        print("pdf(\"%s\")" % (outfile + '.clipping_profile.pdf'), file=ROUT)
                        print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
                        print("clip_count=c(%s)" % ','.join([str(i) for i in clip_count]), file=ROUT)
                        print("nonclip_count= %d - clip_count" % (total_read), file=ROUT)
                        print('plot(read_pos, nonclip_count*100/(clip_count+nonclip_count),col="blue",main="clipping profile",xlab="Position of read",ylab="Non-clipped %",type="b")', file=ROUT)
                        print("dev.off()", file=ROUT)

                
                if PE is True:
                        total_read1 = 0.0
                        total_read2 = 0.0
                        r1_soft_clip_profile = collections.defaultdict(int)
                        r2_soft_clip_profile = collections.defaultdict(int)
                        try:
                                while(1):
                                        aligned_read = next(self.samfile)
                                        if aligned_read.mapq < q_cut: continue
                                        if aligned_read.is_unmapped:continue    #skip unmapped read
                                        if aligned_read.is_qcfail:continue              #skip low quality
                                        if not aligned_read.is_paired: continue
                                        if aligned_read.is_read1:
                                                total_read1 += 1
                                        if aligned_read.is_read2:
                                                total_read2 += 1                                                                                        
                                        cigar_str = bam_cigar.list2longstr(aligned_read.cigar)  # ([(0, 9), (4, 1)] ==> MMMMMMMMMS
                                        if aligned_read.is_reverse:
                                                cigar_str = cigar_str[::-1]                                     
                                        
                                        if type not in cigar_str:       # no clipping
                                                continue
                                        
                                        if aligned_read.is_read1:
                                                for indx,symbl in enumerate(cigar_str):
                                                        if symbl == type:
                                                                r1_soft_clip_profile[indx] += 1.0
                                        if aligned_read.is_read2:
                                                for indx,symbl in  enumerate(cigar_str):
                                                        if symbl == type:
                                                                r2_soft_clip_profile[indx] += 1.0
                        except StopIteration:
                                print("Done", file=sys.stderr)
                        
                        
                        read_pos = list(range(0,len(cigar_str)))
                        r1_clip_count = []
                        r2_clip_count = []
                        
                        print("Totoal read-1 used: %d" % int(total_read1), file=sys.stderr)
                        print("Totoal read-2 used: %d" % int(total_read2), file=sys.stderr)
                        print("Read-1:", file=OUT)
                        for i in read_pos:
                                print(str(i) + '\t' + str(r1_soft_clip_profile[i]) + '\t' + str(total_read1 - r1_soft_clip_profile[i]), file=OUT)               
                                r1_clip_count.append(r1_soft_clip_profile[i])
                                

                        print("Read-2:", file=OUT)
                        for i in read_pos:
                                print(str(i) + '\t' + str(r2_soft_clip_profile[i]) + '\t' + str(total_read2 - r2_soft_clip_profile[i]), file=OUT)
                                r2_clip_count.append(r2_soft_clip_profile[i])
                
                        
                        print("pdf(\"%s\")" % (outfile + '.clipping_profile.R1.pdf'), file=ROUT)
                        print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
                        print("r1_clip_count=c(%s)" % ','.join([str(i) for i in r1_clip_count]), file=ROUT)
                        print("r1_nonclip_count = %d - r1_clip_count" % (total_read1), file=ROUT)
                        print('plot(read_pos, r1_nonclip_count*100/(r1_clip_count + r1_nonclip_count),col="blue",main="clipping profile",xlab="Position of read (read-1)",ylab="Non-clipped %",type="b")', file=ROUT)
                        print("dev.off()", file=ROUT)

                        print("pdf(\"%s\")" % (outfile + '.clipping_profile.R2.pdf'), file=ROUT)
                        print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
                        print("r2_clip_count=c(%s)" % ','.join([str(i) for i in r2_clip_count]), file=ROUT)
                        print("r2_nonclip_count = %d - r2_clip_count" % (total_read2), file=ROUT)
                        print('plot(read_pos, r2_nonclip_count*100/(r2_clip_count + r2_nonclip_count),col="blue",main="clipping profile",xlab="Position of read (read-2)",ylab="Non-clipped %",type="b")', file=ROUT)
                        print("dev.off()", file=ROUT)
                        
        def insertion_profile(self,outfile, q_cut, PE, type="I"):
                '''calculate profile of insertion'''
                
                out_file1 = outfile + ".insertion_profile.xls"
                out_file2 = outfile + ".insertion_profile.r"
                OUT = open(out_file1,'w')
                ROUT = open(out_file2,'w')
                
                print("Position\tInsert_nt\tNon_insert_nt", file=OUT)
                
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)
                
                cigar_str = ""
                
                #single end sequencing
                if PE is False:
                        total_read = 0.0                
                        soft_clip_profile = collections.defaultdict(int)
                        try:
                                while(1):
                                        aligned_read = next(self.samfile)
                                        if aligned_read.mapq < q_cut: continue
                                        if aligned_read.is_unmapped:continue    #skip unmapped read
                                        if aligned_read.is_qcfail:continue              #skip low quality
                                
                                        total_read +=1
                                        cigar_str = bam_cigar.list2longstr(aligned_read.cigar)  # ([(0, 9), (4, 1)] ==> MMMMMMMMMS
                                        
                                        if type not in cigar_str:       # no insertion
                                                continue
                                        if aligned_read.is_reverse:
                                                cigar_str = cigar_str[::-1]             

                                        for indx,symbl in enumerate(cigar_str):
                                                if symbl == type:
                                                        soft_clip_profile[indx] += 1.0                  
                        except StopIteration:
                                print("Done", file=sys.stderr)
                        
                        print("Totoal reads used: %d" % int(total_read), file=sys.stderr)
                        read_pos = list(range(0,len(cigar_str)))
                        clip_count = []
                        for i in read_pos:
                                print(str(i) + '\t' + str(soft_clip_profile[i]) + '\t' + str(total_read - soft_clip_profile[i]), file=OUT)
                                clip_count.append(soft_clip_profile[i])
                        
                        
                        print("pdf(\"%s\")" % (outfile + '.insertion_profile.pdf'), file=ROUT)
                        print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
                        print("insert_count=c(%s)" % ','.join([str(i) for i in clip_count]), file=ROUT)
                        print("noninsert_count= %d - insert_count" % (total_read), file=ROUT)
                        print('plot(read_pos, insert_count*100/(insert_count+noninsert_count),col="blue",main="Insertion profile",xlab="Position of read",ylab="Insertion %",type="b")', file=ROUT)
                        print("dev.off()", file=ROUT)

                
                if PE is True:
                        total_read1 = 0.0
                        total_read2 = 0.0
                        r1_soft_clip_profile = collections.defaultdict(int)
                        r2_soft_clip_profile = collections.defaultdict(int)
                        try:
                                while(1):
                                        aligned_read = next(self.samfile)
                                        if aligned_read.mapq < q_cut: continue
                                        if aligned_read.is_unmapped:continue    #skip unmapped read
                                        if aligned_read.is_qcfail:continue              #skip low quality
                                        if not aligned_read.is_paired: continue
                                        if aligned_read.is_read1:
                                                total_read1 += 1
                                        if aligned_read.is_read2:
                                                total_read2 += 1                                                                                        
                                        cigar_str = bam_cigar.list2longstr(aligned_read.cigar)  # ([(0, 9), (4, 1)] ==> MMMMMMMMMS
                                        if aligned_read.is_reverse:
                                                cigar_str = cigar_str[::-1]                                     
                                        
                                        if type not in cigar_str:       # no clipping
                                                continue
                                        
                                        if aligned_read.is_read1:
                                                for indx,symbl in enumerate(cigar_str):
                                                        if symbl == type:
                                                                r1_soft_clip_profile[indx] += 1.0
                                        if aligned_read.is_read2:
                                                for indx,symbl in  enumerate(cigar_str):
                                                        if symbl == type:
                                                                r2_soft_clip_profile[indx] += 1.0
                        except StopIteration:
                                print("Done", file=sys.stderr)
                        
                        
                        read_pos = list(range(0,len(cigar_str)))
                        r1_clip_count = []
                        r2_clip_count = []
                        
                        print("Totoal read-1 used: %d" % int(total_read1), file=sys.stderr)
                        print("Totoal read-2 used: %d" % int(total_read2), file=sys.stderr)
                        print("Read-1:", file=OUT)
                        for i in read_pos:
                                print(str(i) + '\t' + str(r1_soft_clip_profile[i]) + '\t' + str(total_read1 - r1_soft_clip_profile[i]), file=OUT)               
                                r1_clip_count.append(r1_soft_clip_profile[i])
                                

                        print("Read-2:", file=OUT)
                        for i in read_pos:
                                print(str(i) + '\t' + str(r2_soft_clip_profile[i]) + '\t' + str(total_read2 - r2_soft_clip_profile[i]), file=OUT)
                                r2_clip_count.append(r2_soft_clip_profile[i])
                
                        
                        print("pdf(\"%s\")" % (outfile + '.insertion_profile.R1.pdf'), file=ROUT)
                        print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
                        print("r1_insert_count=c(%s)" % ','.join([str(i) for i in r1_clip_count]), file=ROUT)
                        print("r1_noninsert_count = %d - r1_insert_count" % (total_read1), file=ROUT)
                        print('plot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col="blue",main="Insertion profile",xlab="Position of read (read-1)",ylab="Insertion %",type="b")', file=ROUT)
                        print("dev.off()", file=ROUT)

                        print("pdf(\"%s\")" % (outfile + '.insertion_profile.R2.pdf'), file=ROUT)
                        print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
                        print("r2_insert_count=c(%s)" % ','.join([str(i) for i in r2_clip_count]), file=ROUT)
                        print("r2_noninsert_count = %d - r2_insert_count" % (total_read2), file=ROUT)
                        print('plot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col="blue",main="Insertion profile",xlab="Position of read (read-2)",ylab="Insertion %",type="b")', file=ROUT)
                        print("dev.off()", file=ROUT)
        
        def coverageGeneBody(self,refbed,outfile):
                '''Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divided
                into 100 regsions'''
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                OUT1 = open(outfile + ".geneBodyCoverage_plot.r",'w')
                OUT2 = open(outfile + ".geneBodyCoverage.txt",'w')

                ranges={}
                totalReads=0
                fragment_num=0          #splice reads will counted twice
                rpkm={}
                
                #read SAM 
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                totalReads +=1

                                chrom = self.samfile.getrname(aligned_read.tid).upper()
                                hit_st = aligned_read.pos               
                                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)                   
                                fragment_num += len(exon_blocks)
                        
                                for exon in exon_blocks:
                                        if chrom not in ranges:
                                                ranges[chrom] = Intersecter()
                                        ranges[chrom].add_interval( Interval( exon[1], exon[2] ) )
                except StopIteration:
                        print("Done", file=sys.stderr)

                print("calculating coverage over gene body ...", file=sys.stderr)
                coverage=collections.defaultdict(int)
                flag=0
                for line in open(refbed,'r'):
                        try:
                                if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                                fields = line.split()
                                chrom     = fields[0].upper()
                                tx_start  = int( fields[1] )
                                tx_end    = int( fields[2] )
                                geneName      = fields[3]
                                strand    = fields[5]
                                
                                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        except:
                                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                                continue
                        gene_all_base=[]
                        percentile_base=[]
                        mRNA_len =0
                        flag=0
                        for st,end in zip(exon_starts,exon_ends):
                                gene_all_base.extend(list(range(st+1,end+1)))           #0-based coordinates on genome
                                mRNA_len = len(gene_all_base)
                        if mRNA_len <100:
                                continue
                        if strand == '-':
                                gene_all_base.sort(reverse=True)                        #deal with gene on minus stand
                        else:
                                gene_all_base.sort(reverse=False)
                        percentile_base = mystat.percentile_list (gene_all_base)        #get 101 points from each gene's coordinates
                        
                        for i in range(0,len(percentile_base)):
                                if chrom in ranges:
                                        coverage[i] += len(ranges[chrom].find(percentile_base[i], percentile_base[i]+1))
                x_coord=[]
                y_coord=[]
                print("Total reads: " + str(totalReads), file=OUT2)
                print("Fragment number: " + str(fragment_num), file=OUT2)
                print("percentile\tcount", file=OUT2)
                for i in coverage:
                        x_coord.append(str(i))
                        y_coord.append(str(coverage[i]))
                        print(str(i) + '\t' + str(coverage[i]), file=OUT2)
                print("pdf(\'%s\')" % (outfile + ".geneBodyCoverage.pdf"), file=OUT1)
                print("x=0:100", file=OUT1)
                print("y=c(" + ','.join(y_coord) + ')', file=OUT1)
                print("plot(x,y,xlab=\"percentile of gene body (5'->3')\",ylab='read number',type='s')", file=OUT1)
                print("dev.off()", file=OUT1)
                                        
        def mRNA_inner_distance(self,outfile,refbed,low_bound=0,up_bound=1000,step=10,sample_size=1000000, q_cut=30):
                '''estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length'''
                
                out_file1 = outfile + ".inner_distance.txt"     
                out_file2 = outfile + ".inner_distance_freq.txt"
                out_file3 = outfile + ".inner_distance_plot.r"
                
                FO=open(out_file1,'w')
                FQ=open(out_file2,'w')
                RS=open(out_file3,'w')
                
                fchrom="chr100"         #this is the fake chromosome
                ranges={}
                ranges[fchrom]=Intersecter()
                
                window_left_bound = list(range(low_bound,up_bound,step))
                frag_size=0
                
                inner_distance_bitsets=BinnedBitSet()
                tmp = BinnedBitSet()
                tmp.set_range(0,0)
                pair_num=0
                sizes=[]
                counts=[]
                count=0
                
                print("Get exon regions from " + refbed + " ...", file=sys.stderr)
                bed_obj = BED.ParseBED(refbed)
                ref_exons = []
                
                for exn in bed_obj.getExon():
                        ref_exons.append([exn[0].upper(), exn[1], exn[2]])
                exon_bitsets = binned_bitsets_from_list(ref_exons)
                
                transcript_ranges = {}
                for i_chr, i_st, i_end, i_strand, i_name in bed_obj.getTranscriptRanges():
                        i_chr = i_chr.upper()
                        if i_chr not in transcript_ranges:
                                transcript_ranges[i_chr] = Intersecter()
                        else:
                                transcript_ranges[i_chr].add_interval(Interval(i_st, i_end, value=i_name))
                
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                if pair_num >= sample_size:
                                        break
                                splice_intron_size=0
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                if not aligned_read.is_paired: continue         #skip single map read
                                if aligned_read.mate_is_unmapped:continue       #
                                if aligned_read.mapq < q_cut:continue
                                
                                read1_len = aligned_read.qlen
                                read1_start = aligned_read.pos
                                read2_start = aligned_read.mpos         #0-based, not included
                                if read2_start < read1_start:
                                        continue                                                                #because BAM file is sorted, mate_read is already processed if its coordinate is smaller
                                if  read2_start == read1_start and aligned_read.is_read1:
                                        inner_distance = 0
                                        continue
                                
                                pair_num +=1

                                # check if reads were mapped to diff chromsomes
                                R_read1_ref = self.samfile.getrname(aligned_read.tid)
                                R_read2_ref = self.samfile.getrname(aligned_read.rnext)
                                if R_read1_ref != R_read2_ref:
                                        FO.write(aligned_read.qname + '\t' + 'NA' + '\tsameChrom=No\n') #reads mapped to different chromosomes
                                        continue
                                
                                
                                chrom = self.samfile.getrname(aligned_read.tid).upper()
                                intron_blocks = bam_cigar.fetch_intron(chrom, read1_start, aligned_read.cigar)                          
                                for intron in intron_blocks:
                                        splice_intron_size += intron[2] - intron[1]
                                read1_end = read1_start + read1_len + splice_intron_size                
                                
                                if read2_start >= read1_end:
                                        inner_distance = read2_start - read1_end
                                else:
                                        exon_positions = []
                                        exon_blocks = bam_cigar.fetch_exon(chrom, read1_start,aligned_read.cigar)
                                        for ex in exon_blocks:
                                                for i in range(ex[1]+1,ex[2]+1):
                                                        exon_positions.append(i)
                                        inner_distance = -len([i for i in exon_positions if i > read2_start and i <= read1_end])
                                #print aligned_read.qname,read1_end, read2_start
                                
                                read1_gene_names = set()        #read1_end
                                try:
                                        for gene in transcript_ranges[chrom].find(read1_end-1, read1_end):      #gene: Interval(0, 10, value=a)
                                                read1_gene_names.add(gene.value)
                                except:
                                        pass
                                
                                read2_gene_names = set()        #read2_start
                                try:
                                        for gene in transcript_ranges[chrom].find(read2_start, read2_start +1): #gene: Interval(0, 10, value=a)
                                                read2_gene_names.add(gene.value)
                                except:
                                        pass
                                        
                                if len(read1_gene_names.intersection(read2_gene_names)) == 0:   # no common gene
                                        FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tsameTranscript=No,dist=genomic\n')                #reads mapped to different gene
                                        ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )             
                                        continue                        

                                if inner_distance > 0: 
                                        if chrom in exon_bitsets:
                                                size =0 
                                                inner_distance_bitsets.set_range(read1_end, read2_start-read1_end)
                                                inner_distance_bitsets.iand(exon_bitsets[chrom])
                                                end=0
                                                while 1:
                                                        start = inner_distance_bitsets.next_set( end )
                                                        if start == inner_distance_bitsets.size: break
                                                        end = inner_distance_bitsets.next_clear( start )
                                                        size += (end - start)
                                                inner_distance_bitsets.iand(tmp)                                                                                        #clear BinnedBitSet
                                                
                                                if size == inner_distance:
                                                        FO.write(aligned_read.qname + '\t' + str(size) + '\tsameTranscript=Yes,sameExon=Yes,dist=mRNA\n')
                                                        ranges[fchrom].add_interval( Interval( size-1, size ) )
                                                elif size > 0 and size < inner_distance:
                                                        FO.write(aligned_read.qname + '\t' + str(size) + '\tsameTranscript=Yes,sameExon=No,dist=mRNA\n')
                                                        ranges[fchrom].add_interval( Interval( size-1, size ) ) 
                                                elif size <= 0:
                                                        FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tsameTranscript=Yes,nonExonic=Yes,dist=genomic\n')
                                                        ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )             
                                        else:
                                                FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tunknownChromosome,dist=genomic')
                                                ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
                                else:
                                        FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\treadPairOverlap\n')
                                        ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
                                
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                print("Total read pairs  used " + str(pair_num), file=sys.stderr)
                if pair_num==0:
                        print("Cannot find paired reads", file=sys.stderr)
                        sys.exit(0)
                
                for st in window_left_bound:
                        sizes.append(str(st + step/2))
                        count = str(len(ranges[fchrom].find(st,st + step)))
                        counts.append(count)
                        print(str(st) + '\t' + str(st+step) +'\t' + count, file=FQ)             
                print("out_file = \'%s\'" % outfile, file=RS)
                print("pdf(\'%s\')" % (outfile + ".inner_distance_plot.pdf"), file=RS)
                #print >>RS, "par(mfrow=c(2,1),cex.main=0.8,cex.lab=0.8,cex.axis=0.8,mar=c(4,4,4,1))"
                #print >>RS, 'pie(c(%d,%d,%d),col=rainbow(3),cex=0.5,radius=1,main="Total %d fragments",labels=c("fraSize <= %d\\n(%4.2f%%)","fragSize > %d\\n(%4.2f%%)","%d < fragSize <= %d\\n(%4.2f%%)"), density=rep(80,80,80),angle=c(90,140,170))' % (ultra_low, ultra_high, pair_num -ultra_low -ultra_high, pair_num, low_bound, ultra_low*100/pair_num, up_bound, ultra_high*100/pair_num, low_bound, up_bound, 100-ultra_low*100/pair_num - ultra_high*100/pair_num)
                print('fragsize=rep(c(' + ','.join(sizes) + '),' + 'times=c(' + ','.join(counts) + '))', file=RS)
                print('frag_sd = sd(fragsize)', file=RS)
                print('frag_mean = mean(fragsize)', file=RS)
                print('frag_median = median(fragsize)', file=RS)
                print('write(x=c("Name","Mean","Median","sd"), sep="\t", file=stdout(),ncolumns=4)', file=RS)
                print('write(c(out_file,frag_mean,frag_median,frag_sd),sep="\t", file=stdout(),ncolumns=4)', file=RS)
                print('hist(fragsize,probability=T,breaks=%d,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(window_left_bound), file=RS)
                print("lines(density(fragsize,bw=%d),col='red')" % (2*step), file=RS)
                print("dev.off()", file=RS)
                FO.close()
                FQ.close()
                RS.close()
                #self.f.seek(0)
        def annotate_junction(self,refgene,outfile,min_intron=50, q_cut=30):
                '''Annotate splicing junctions in BAM or SAM file. Note that a (long) read might have multiple splicing
                events  (splice multiple times), and the same splicing events can be consolidated into a single
                junction'''
                
                out_file = outfile + ".junction.xls"
                out_file2 = outfile + ".junction_plot.r"
                if refgene is None:
                        print("You must provide reference gene model in bed format.", file=sys.stderr)
                        sys.exit(1)
                OUT = open(out_file,'w')
                ROUT = open(out_file2,'w')
                
                #reading reference gene model
                refIntronStarts=collections.defaultdict(dict)
                refIntronEnds=collections.defaultdict(dict)     
                total_junc =0
                novel35_junc =0
                novel3or5_junc =0
                known_junc =0
                splicing_events=collections.defaultdict(int)    
                
                print("Reading reference bed file: ",refgene, " ... ", end=' ', file=sys.stderr)
                for line in open(refgene,'r'):
                        if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                        fields = line.split()
                        if(len(fields)<12):
                                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                                continue
                        chrom     = fields[0].upper()
                        tx_start = int( fields[1] )
                        tx_end   = int( fields[2] )
                        if int(fields[9] ==1):
                                continue        
                        
                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        intron_start = exon_ends[:-1]
                        intron_end=exon_starts[1:]
                        for i_st,i_end in zip (intron_start, intron_end):
                                refIntronStarts[chrom][i_st] =i_st
                                refIntronEnds[chrom][i_end] =i_end                      
                print("Done", file=sys.stderr)
                
                #reading input SAM file
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)

                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                if aligned_read.mapq < q_cut:continue

                                chrom = self.samfile.getrname(aligned_read.tid).upper()
                                hit_st = aligned_read.pos
                                intron_blocks = bam_cigar.fetch_intron(chrom, hit_st, aligned_read.cigar)                       
                                if len(intron_blocks)==0:
                                        continue
                                for intrn in intron_blocks:
                                        total_junc +=1
                                        if intrn[2] - intrn[1] < min_intron:continue
                                        splicing_events[intrn[0] + ":" + str(intrn[1]) + ":" + str(intrn[2])] += 1
                                        if (intrn[1] in refIntronStarts[chrom] and intrn[2] in refIntronEnds[chrom]):
                                                known_junc +=1                                                                                                                                          #known both
                                        elif (intrn[1] not in refIntronStarts[chrom] and intrn[2] not in refIntronEnds[chrom]):
                                                novel35_junc +=1                                                                                                                                
                                        else:
                                                novel3or5_junc +=1
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                print("total = " + str(total_junc))
                if total_junc == 0:
                        print("No splice junction found.", file=sys.stderr)
                        sys.exit()
                #self.f.seek(0)
                
                print('pdf(\"%s\")' % (outfile + ".splice_events.pdf"), file=ROUT)
                print("events=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc)])+ ')', file=ROUT)
                print('pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing events",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc)), file=ROUT)
                print("dev.off()", file=ROUT)
                
                print("\n===================================================================", file=sys.stderr)
                print("Total splicing  Events:\t" + str(total_junc), file=sys.stderr)
                print("Known Splicing Events:\t" + str(known_junc), file=sys.stderr)
                print("Partial Novel Splicing Events:\t" + str(novel3or5_junc), file=sys.stderr)
                print("Novel Splicing Events:\t" + str(novel35_junc), file=sys.stderr)
                
                #reset variables
                total_junc =0
                novel35_junc =0
                novel3or5_junc =0
                known_junc =0
                
                print("chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation", file=OUT)
                for i in splicing_events:
                        total_junc += 1
                        (chrom, i_st, i_end) = i.split(":")
                        print('\t'.join([chrom.replace("CHR","chr"),i_st,i_end]) + '\t' + str(splicing_events[i]) + '\t', end=' ', file=OUT)
                        i_st = int(i_st)
                        i_end = int(i_end)
                        if (i_st in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]):
                                print("annotated", file=OUT)
                                known_junc +=1
                        elif (i_st not in refIntronStarts[chrom] and i_end not in refIntronEnds[chrom]):
                                print('complete_novel', file=OUT)
                                novel35_junc +=1
                        else:
                                print('partial_novel', file=OUT)
                                novel3or5_junc +=1
                
                if total_junc ==0:
                        print("No splice read found", file=sys.stderr)
                        sys.exit(1)
                print("\nTotal splicing  Junctions:\t" + str(total_junc), file=sys.stderr)
                print("Known Splicing Junctions:\t" + str(known_junc), file=sys.stderr)
                print("Partial Novel Splicing Junctions:\t" + str(novel3or5_junc), file=sys.stderr)
                print("Novel Splicing Junctions:\t" + str(novel35_junc), file=sys.stderr)
                print("\n===================================================================", file=sys.stderr)
                
                print('pdf(\"%s\")' % (outfile + ".splice_junction.pdf"), file=ROUT)
                print("junction=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc,)])+ ')', file=ROUT)
                print('pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing junctions",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc)), file=ROUT)
                print("dev.off()", file=ROUT)
                #print >>ROUT, "mat=matrix(c(events,junction),byrow=T,ncol=3)"
                #print >>ROUT, 'barplot(mat,beside=T,ylim=c(0,100),names=c("known","partial\nnovel","complete\nnovel"),legend.text=c("splicing events","splicing junction"),ylab="Percent")'
        def junction_freq(self, chrom, st, end, known_junctions, q_cut=30):
                '''
                return number of splicing reads for each known junction
                '''
                
                #reading input SAM file
                #if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
                #else:print >>sys.stderr, "Load SAM file ... ",
                
                junc_freq = collections.defaultdict(int)
                try:
                        alignedReads = self.samfile.fetch(chrom,st,end)
                except:
                        return junc_freq
                for aligned_read in alignedReads:
                        if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                        if aligned_read.is_duplicate:continue           #skip duplicate read
                        if aligned_read.is_secondary:continue           #skip non primary hit
                        if aligned_read.is_unmapped:continue            #skip unmap read
                        if aligned_read.mapq < q_cut:continue

                        intron_blocks = bam_cigar.fetch_intron(chrom, aligned_read.pos, aligned_read.cigar)                     
                        if len(intron_blocks)==0:
                                continue
                        for intrn in intron_blocks:
                                tmp = chrom + ":" + str(intrn[1]) + '-' + str(intrn[2])
                                if tmp in known_junctions:
                                        junc_freq[tmp] += 1
                                else:
                                        continue
                for k in known_junctions:
                        if k not in list(junc_freq.keys()):
                                junc_freq[k] = 0
                        elif junc_freq[k] < 2: junc_freq[k] = 0
                return junc_freq
        def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1, q_cut=30):
                '''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction'''
                
                out_file = outfile + ".junctionSaturation_plot.r"
                if refgene is None:
                        print("You must provide reference gene model in bed format.", file=sys.stderr)
                        sys.exit(1)
                
                OUT = open(out_file,'w')


                #reading reference gene 
                knownSpliceSites= set()
                chrom_list=set()
                print("reading reference bed file: ",refgene, " ... ", end=' ', file=sys.stderr)
                for line in open(refgene,'r'):
                        if line.startswith(('#','track','browser')):continue  
                        fields = line.split()
                        if(len(fields)<12):
                                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                                continue
                        chrom     = fields[0].upper()
                        chrom_list.add(chrom)
                        tx_start = int( fields[1] )
                        tx_end   = int( fields[2] )
                        if int(fields[9] ==1):
                                continue        
                        
                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                        intron_start = exon_ends[:-1]
                        intron_end=exon_starts[1:]
                        for st,end in zip (intron_start, intron_end):
                                knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))
                print("Done! Total "+str(len(knownSpliceSites)) + " known splicing junctions.", file=sys.stderr)


                #read SAM file
                samSpliceSites=[]
                intron_start=[]
                intron_end=[]
                uniqSpliceSites=collections.defaultdict(int)

                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)
                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                try:
                                        chrom = self.samfile.getrname(aligned_read.tid).upper()
                                except:
                                        continue
                                if chrom not in chrom_list:
                                        continue                                
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                if aligned_read.mapq < q_cut: continue
                                
                                hit_st = aligned_read.pos
                                intron_blocks = bam_cigar.fetch_intron(chrom, hit_st, aligned_read.cigar)                       
                                if len(intron_blocks)==0:
                                        continue
                                for intrn in intron_blocks:
                                        if intrn[2] - intrn[1] < min_intron:continue
                                        samSpliceSites.append(intrn[0] + ":" + str(intrn[1]) + "-" + str(intrn[2]))
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                print("shuffling alignments ...", end=' ', file=sys.stderr)
                random.shuffle(samSpliceSites)
                print("Done", file=sys.stderr)
                                
                #resampling
                SR_num = len(samSpliceSites)
                sample_size=0
                all_junctionNum = 0     
                known_junc=[]
                all_junc=[]
                unknown_junc=[]
                #=========================sampling uniquely mapped reads from population
                tmp=list(range(sample_start,sample_end,sample_step))
                tmp.append(100)
                for pertl in tmp:       #[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
                        knownSpliceSites_num = 0
                        index_st = int(SR_num * ((pertl - sample_step)/100.0))
                        index_end = int(SR_num * (pertl/100.0))
                        if index_st < 0: index_st = 0
                        sample_size += index_end -index_st
                        
                        print("sampling " + str(pertl) +"% (" + str(sample_size) + ") splicing reads.", end=' ', file=sys.stderr)
                        
                        #all splice juntion
                        for i in range(index_st, index_end):
                                uniqSpliceSites[samSpliceSites[i]] +=1  
                        all_junctionNum = len(list(uniqSpliceSites.keys()))
                        all_junc.append(str(all_junctionNum))
                        print(str(all_junctionNum) + " splicing junctions.", end=' ', file=sys.stderr)
                        
                        #known splice junction
                        known_junctionNum = 0
                        for sj in uniqSpliceSites:
                                if sj in knownSpliceSites and uniqSpliceSites[sj] >= recur:
                                        known_junctionNum +=1
                        print(str(known_junctionNum) + " known splicing junctions.", end=' ', file=sys.stderr)
                        known_junc.append(str(known_junctionNum))
                        
                        #unknown splice junction
                        unknown_junctionNum = 0
                        for sj in uniqSpliceSites:
                                if sj not in knownSpliceSites:
                                        unknown_junctionNum +=1
                        unknown_junc.append(str(unknown_junctionNum))
                        print(str(unknown_junctionNum) + " novel splicing junctions.", file=sys.stderr)
                        
                #for j in uniq_SJ:
                        #print >>OUT, j + "\t" + str(uniq_SJ[j])
                print("pdf(\'%s\')" % (outfile + '.junctionSaturation_plot.pdf'), file=OUT)
                print("x=c(" + ','.join([str(i) for i in tmp]) + ')', file=OUT)
                print("y=c(" + ','.join(known_junc) + ')', file=OUT)
                print("z=c(" + ','.join(all_junc) + ')', file=OUT)
                print("w=c(" + ','.join(unknown_junc) + ')', file=OUT)
                print("m=max(%d,%d,%d)" % (int(int(known_junc[-1])/1000), int(int(all_junc[-1])/1000),int(int(unknown_junc[-1])/1000)), file=OUT)
                print("n=min(%d,%d,%d)" % (int(int(known_junc[0])/1000), int(int(all_junc[0])/1000),int(int(unknown_junc[0])/1000)), file=OUT)
                print("plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))", file=OUT)
                print("points(x,y/1000,type='o',col='red')", file=OUT)
                print("points(x,w/1000,type='o',col='green')", file=OUT)
                print('legend(5,%d, legend=c("All junctions","known junctions", "novel junctions"),col=c("blue","red","green"),lwd=1,pch=1)' % int(int(all_junc[-1])/1000), file=OUT)
                print("dev.off()", file=OUT)
        def saturation_RPKM(self,refbed,outfile,sample_start=5,sample_step=5,sample_end=100,skip_multi=True, strand_rule=None, q_cut=30):
                '''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                rpkm_file = outfile + ".eRPKM.xls"
                raw_file = outfile + ".rawCount.xls"
                
                RPKM_OUT = open(rpkm_file,'w')
                RAW_OUT = open(raw_file ,'w')
                
                ranges={}
                totalReads=0
                cUR_num = 0     #number of fragements
                cUR_plus = 0
                cUR_minus = 0
                block_list_plus = []    #non-spliced read AS IS, splicing reads were counted multiple times
                block_list_minus = []
                block_list = []
                strandRule = {}
                                
                if strand_rule is None:                                                                                                 # Not strand-specific
                        pass                                                                                                                            
                elif len(strand_rule.split(',')) ==4:                                                                   #PairEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
                elif len(strand_rule.split(',')) ==2:                                                                   #singeEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]]=i[1]
                else:
                        print("Unknown value of: 'strand_rule' " +  strand_rule, file=sys.stderr)
                        sys.exit(1)     


                #read SAM or BAM
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)
                try:
                        while(1):
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                
                                if skip_multi:
                                        if aligned_read.mapq < q_cut:
                                                continue                        
                                chrom = self.samfile.getrname(aligned_read.tid).upper()
                                
                                #determine read_id and read_strand
                                if aligned_read.is_paired:                                              #pair end
                                        if aligned_read.is_read1:read_id = '1'
                                        if aligned_read.is_read2:read_id = '2'
                                else:read_id = ''                                                               #single end
                        
                                if aligned_read.is_reverse:map_strand = '-'
                                else:map_strand = '+'                           
                                strand_key = read_id + map_strand                               #used to determine if a read should assign to gene(+) or gene(-)

                                hit_st = aligned_read.pos
                                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)   
                                cUR_num += len(exon_blocks)
                                
                                #strand specific
                                if strand_rule is not None:
                                        if strandRule[strand_key] == '+': cUR_plus += len(exon_blocks)
                                        if strandRule[strand_key] == '-': cUR_minus += len(exon_blocks)
                                        for exn in exon_blocks:
                                                if strandRule[strand_key] == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
                                                if strandRule[strand_key] == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
                                #Not strand specific
                                else:                   
                                        for exn in exon_blocks:
                                                block_list.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
                except StopIteration:
                        print("Done", file=sys.stderr)
                
                
                print("shuffling alignments ...", end=' ', file=sys.stderr)
                random.shuffle(block_list_plus)
                random.shuffle(block_list_minus)
                random.shuffle(block_list)
                print("Done", file=sys.stderr)
                
                
                ranges_plus={}
                ranges_minus={}
                ranges={}
                sample_size=0
                RPKM_table=collections.defaultdict(list)
                rawCount_table=collections.defaultdict(list)
                RPKM_head=['#chr','start','end','name','score','strand']

                tmp=list(range(sample_start,sample_end,sample_step))
                tmp.append(100)
                #=========================sampling uniquely mapped reads from population
                for pertl in tmp:       #[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
                        percent_st = (pertl-sample_step)/100.0
                        percent_end = pertl/100.0
                        if percent_st < 0: percent_st = 0
                        sample_size = cUR_num * percent_end
                        RPKM_head.append(str(pertl) + '%')
                        
                        if strand_rule is not None:
                                print("sampling " + str(pertl) +"% (" + str(int(cUR_plus * percent_end)) + ") forward strand fragments ...", file=sys.stderr)
                                for i in block_list_plus[int(cUR_plus*percent_st):int(cUR_plus*percent_end)]:
                                        (chr,coord) = i.split(':')
                                        if chr not in ranges_plus:ranges_plus[chr] = Intersecter()
                                        ranges_plus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                                           
                                
                                print("sampling " + str(pertl) +"% (" + str(int(cUR_minus * percent_end)) + ") reverse strand fragments ...", file=sys.stderr)                  
                                for i in block_list_minus[int(cUR_minus*percent_st):int(cUR_minus*percent_end)]:
                                        (chr,coord) = i.split(':')
                                        if chr not in ranges_minus:ranges_minus[chr] = Intersecter()                            
                                        ranges_minus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                          
                        
                        else:
                                print("sampling " + str(pertl) +"% (" + str(int(sample_size)) + ") fragments ...", file=sys.stderr)
                                for i in block_list[int(cUR_num*percent_st):int(cUR_num*percent_end)]:
                                        (chr,coord) = i.split(':')
                                        if chr not in ranges:ranges[chr] = Intersecter()                                                
                                        ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                                                                                                

                        #========================= calculating RPKM based on sub-population
                        print("assign reads to transcripts in " + refbed + ' ...', file=sys.stderr)
                        for line in open(refbed,'r'):
                                try:
                                        if line.startswith(('#','track','browser')):continue  
                        # Parse fields from gene tabls
                                        fields = line.split()
                                        chrom     = fields[0].upper()
                                        tx_start  = int( fields[1] )
                                        tx_end    = int( fields[2] )
                                        geneName      = fields[3]
                                        strand    = fields[5]
                                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                        exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                        key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                                except:
                                        print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                                        continue
                                mRNA_count=0    #we need to initializ it to 0 for each gene
                                mRNA_len=sum(exon_sizes)
                                for st,end in zip(exon_starts,exon_ends):
                                        #if chrom in ranges:
                                        if strand_rule is not None:
                                                if (strand == '+') and (chrom in ranges_plus): mRNA_count += len(ranges_plus[chrom].find(st,end))       
                                                if (strand == '-') and (chrom in ranges_minus): mRNA_count += len(ranges_minus[chrom].find(st,end))
                                        else:
                                                if chrom in ranges:
                                                        mRNA_count += len(ranges[chrom].find(st,end))
                                if mRNA_len ==0:
                                        print(geneName + " has 0 nucleotides. Exit!", file=sys.stderr)
                                        sys.exit(1)
                                if sample_size == 0:
                                        print("Too few reads to sample. Exit!", file=sys.stderr)
                                        sys.exit(1)
                                mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * sample_size)
                                RPKM_table[key].append(str(mRNA_RPKM))
                                rawCount_table[key].append(str(mRNA_count))
                        print("", file=sys.stderr)

                #self.f.seek(0)
                print('\t'.join(RPKM_head), file=RPKM_OUT)
                print('\t'.join(RPKM_head), file=RAW_OUT)
                for key in RPKM_table:
                        print(key + '\t', end=' ', file=RPKM_OUT)
                        print('\t'.join(RPKM_table[key]), file=RPKM_OUT)
                        print(key + '\t', end=' ', file=RAW_OUT)
                        print('\t'.join(rawCount_table[key]), file=RAW_OUT)             
        def shuffle_RPKM(self,refbed,outfile,sample_percentage=0.5,shuffle_times=50,skip_multi=True, strand_rule=None):
                '''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
                
                if refbed is None:
                        print("You must specify a bed file representing gene model\n", file=sys.stderr)
                        exit(0)
                rpkm_file = outfile + ".eRPKM.xls"
                raw_file = outfile + ".rawCount.xls"
                
                RPKM_OUT = open(rpkm_file,'w')
                RAW_OUT = open(raw_file ,'w')
                
                ranges={}
                totalReads=0
                cUR_num = 0     #number of fragements
                cUR_plus = 0
                cUR_minus = 0
                block_list_plus = []    #non-spliced read AS IS, splicing reads were counted multiple times
                block_list_minus = []
                block_list = []
                strandRule = {}
                                
                if strand_rule is None:                                                                                                 # Not strand-specific
                        pass                                                                                                                            
                elif len(strand_rule.split(',')) ==4:                                                                   #PairEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
                elif len(strand_rule.split(',')) ==2:                                                                   #singeEnd, strand-specific
                        for i in strand_rule.split(','):strandRule[i[0]]=i[1]
                else:
                        print("Unknown value of: 'strand_rule' " +  strand_rule, file=sys.stderr)
                        sys.exit(1)     


                #read SAM or BAM
                if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
                else:print("Load SAM file ... ", end=' ', file=sys.stderr)
                try:
                        while(1):
                                flag=0
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                
                                if skip_multi:
                                        if len(aligned_read.tags)>0:            #( ("NM", 1),("RG", "L1") )
                                                for i in aligned_read.tags:
                                                        if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
                                                                flag=1                                          #multiple hit read
                                                                break
                                if flag==1:continue                                             #skip multiple map read         
                                
                                chrom = self.samfile.getrname(aligned_read.tid).upper()
                                
                                #determine read_id and read_strand
                                if aligned_read.is_paired:                                              #pair end
                                        if aligned_read.is_read1:read_id = '1'
                                        if aligned_read.is_read2:read_id = '2'
                                else:read_id = ''                                                               #single end
                        
                                if aligned_read.is_reverse:map_strand = '-'
                                else:map_strand = '+'                           
                                strand_key = read_id + map_strand                               #used to determine if a read should assign to gene(+) or gene(-)

                                hit_st = aligned_read.pos
                                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)   
                                cUR_num += len(exon_blocks)
                                
                                #strand specific
                                if strand_rule is not None:
                                        if strandRule[strand_key] == '+': cUR_plus += len(exon_blocks)
                                        if strandRule[strand_key] == '-': cUR_minus += len(exon_blocks)
                                        for exn in exon_blocks:
                                                if strandRule[strand_key] == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
                                                if strandRule[strand_key] == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
                                #Not strand specific
                                else:                   
                                        for exn in exon_blocks:
                                                block_list.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
                except StopIteration:
                        print("Done", file=sys.stderr)
                                
                RPKM_table=collections.defaultdict(list)
                rawCount_table=collections.defaultdict(list)
                RPKM_head=['#chr','start','end','name','score','strand']
                
                iter_times=0
                #=========================sampling uniquely mapped reads from population
                for x in range(0,shuffle_times+1):
                        print("Shuffle " + str(iter_times) + " times", file=sys.stderr)
                        iter_times += 1
                        if iter_times == shuffle_times:
                                sample_percent = 1
                        else:
                                sample_percent = sample_percentage
                        ranges_plus={}
                        ranges_minus={}
                        ranges={}
                        if strand_rule is not None:
                                for i in random.sample(block_list_plus, int(cUR_plus * sample_percent)):
                                        (chr,coord) = i.split(':')
                                        if chr not in ranges_plus:ranges_plus[chr] = Intersecter()
                                        ranges_plus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                                           
                                
                                for i in random.sample(block_list_minus, int(cUR_minus * sample_percent)):
                                        (chr,coord) = i.split(':')
                                        if chr not in ranges_minus:ranges_minus[chr] = Intersecter()                            
                                        ranges_minus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                          
                        
                        else:
                                for i in random.sample(block_list,int(cUR_num * sample_percent)):
                                        (chr,coord) = i.split(':')
                                        if chr not in ranges:ranges[chr] = Intersecter()                                                
                                        ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                                                                                                

                        #========================= calculating RPKM based on sub-population
                        print("assign reads to transcripts in " + refbed + ' ...', file=sys.stderr)
                        for line in open(refbed,'r'):
                                try:
                                        if line.startswith(('#','track','browser')):continue  
                        # Parse fields from gene tabls
                                        fields = line.split()
                                        chrom     = fields[0].upper()
                                        tx_start  = int( fields[1] )
                                        tx_end    = int( fields[2] )
                                        geneName      = fields[3]
                                        strand    = fields[5]
                                        exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                                        exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                                        exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                                        exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                                        exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                                        key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
                                except:
                                        print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                                        continue
                                mRNA_count=0    #we need to initializ it to 0 for each gene
                                mRNA_len=sum(exon_sizes)
                                for st,end in zip(exon_starts,exon_ends):
                                        #if chrom in ranges:
                                        if strand_rule is not None:
                                                if (strand == '+') and (chrom in ranges_plus): mRNA_count += len(ranges_plus[chrom].find(st,end))       
                                                if (strand == '-') and (chrom in ranges_minus): mRNA_count += len(ranges_minus[chrom].find(st,end))
                                        else:
                                                if chrom in ranges:
                                                        mRNA_count += len(ranges[chrom].find(st,end))
                                if mRNA_len ==0:
                                        print(geneName + " has 0 nucleotides. Exit!", file=sys.stderr)
                                        sys.exit(1)
                                if cUR_num * sample_percentage == 0:
                                        print("Too few reads to sample. Exit!", file=sys.stderr)
                                        sys.exit(1)
                                mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * (cUR_num * sample_percentage))
                                RPKM_table[key].append(str(mRNA_RPKM))
                                rawCount_table[key].append(str(mRNA_count))
                        print("", file=sys.stderr)

                #self.f.seek(0)
                print('\t'.join(RPKM_head), file=RPKM_OUT)
                print('\t'.join(RPKM_head), file=RAW_OUT)
                for key in RPKM_table:
                        print(key + '\t', end=' ', file=RPKM_OUT)
                        print('\t'.join(RPKM_table[key]), file=RPKM_OUT)
                        print(key + '\t', end=' ', file=RAW_OUT)
                        print('\t'.join(rawCount_table[key]), file=RAW_OUT)             
        def fetchAlignments(self,chr,st,end):
                '''fetch alignment from sorted BAM file based on chr, st, end
                Note: BAM file must be indexed'''
                try:
                        a=self.samfile.fetch(chr,st,end)
                        return a
                except:
                        return None
        def mismatchProfile(self,read_length,read_num, outfile, q_cut=30):
                '''
                Calculate mismatch profile. Note that the "MD" tag must exist.
                '''
                
                DOUT = open(outfile + '.mismatch_profile.xls','w')
                ROUT = open(outfile + '.mismatch_profile.r','w')
                
                #reading input SAM file
                if self.bam_format:print("Process BAM file ... ", end=' ', file=sys.stderr)
                else:print("Process SAM file ... ", end=' ', file=sys.stderr)
                
                MD_pat =re.compile(r'(\d+)([A-Z]+)')
                number_base = re.compile(r'([0-9]+)([A-Z]+)', re.I)
                
                count = 0
                data = collections.defaultdict(dict)    # data[read_coord][genotype] = geno_type_number
                try:
                        while(1):
                                if count >= read_num:
                                        print("Total reads used: " + str(count), file=sys.stderr)
                                        break
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                if aligned_read.mapq < q_cut: continue
                                if aligned_read.is_reverse:
                                        strand = '-'
                                else:
                                        strand = '+'
                                
                                #Skip if there is no mismatch, or there is deletion
                                tags = aligned_read.tags
                                skip = False
                                for tag in  tags:
                                        if tag[0] == "NM" and tag[1] ==0:
                                                skip = True             #skip reads with no mismatches
                                        if (tag[0] == "MD") and ('^' in tag[1]):
                                                skip = True             # skip as there is deletion from the reference
                                if skip is True: continue
                                
                                # skip partially mapped read
                                read_seq = aligned_read.seq
                                if len(read_seq) != read_length:
                                        continue
                                if 'N' in read_seq:
                                        continue
                                
                                matched_portion_size = 0
                                for op,value in aligned_read.cigar:
                                        if op == 0: matched_portion_size += value
                                if matched_portion_size != read_length:
                                        continue
                                
                                count += 1
                                if strand == '+':
                                        for tag in tags:
                                                if tag[0] == "MD":
                                                        a = MD_pat.findall(tag[1])      #tag[1] = "5G19T75"; a = [('5', 'G'), ('19', 'T')]
                                                        read_coord = 0
                                                        for (match_number, ref_base) in a:
                                                                read_coord += int(match_number)
                                                                read_base = read_seq[read_coord]
                                                                if read_base == ref_base: continue
                                                                genotype = ref_base + '2' + read_base
                                                                if genotype not in data[read_coord]:
                                                                        data[read_coord][genotype] = 1
                                                                else:
                                                                        data[read_coord][genotype] += 1
                                                                read_coord += 1                                 
                                if strand == '-':
                                        for tag in tags:
                                                if tag[0] == "MD":
                                                        a = MD_pat.findall(tag[1])      #tag[1] = "5G19T75"; a = [('5', 'G'), ('19', 'T')]
                                                        read_coord = 0
                                                        for (match_number, ref_base) in a:
                                                                read_coord += int(match_number)
                                                                read_base = read_seq[read_coord]
                                                                if read_base == ref_base: continue
                                                                genotype = ref_base + '2' + read_base
                                                                if genotype not in data[read_length - read_coord -1]:
                                                                        data[read_length - read_coord -1][genotype] = 1
                                                                else:
                                                                        data[read_length - read_coord -1][genotype] += 1
                                                                read_coord += 1
                                                                if read_base == ref_base: print(aligned_read)                           
                
                except StopIteration:
                        print("Total reads used: " + str(count), file=DOUT)
                print('\n')
                
                if len(data) == 0:
                        print("No mismatches found", file=sys.stderr)
                        sys.exit()
                # write data out
                all_genotypes = ['A2C','A2G','A2T','C2A','C2G','C2T','G2A','G2C','G2T','T2A','T2C','T2G']
                print("read_pos\tsum\t" + '\t'.join(all_genotypes), file=DOUT)
                for indx in sorted(data):
                        tmp = [indx, sum(data[indx].values())]  #read position and sum of mismatches
                        for i in all_genotypes:
                                if i in data[indx]:
                                        tmp.append(data[indx][i])
                                else:
                                        tmp.append(0)
                        print('\t'.join([str(i) for i in tmp]), file=DOUT)
                
                DOUT.close()
                
                # write Rscript
                r_data = collections.defaultdict(list)
                for gt in all_genotypes:
                        for indx in sorted(data):
                                if gt in data[indx]:
                                        r_data[gt].append(data[indx][gt])
                                else:
                                        r_data[gt].append(0)
                for k in sorted(r_data):
                        print('%s=c(%s)' % (k, ','.join([str(i) for i in r_data[k]])), file=ROUT)
                
                #print >>ROUT, 'color_code = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")'
                print('color_code = c("green","powderblue","lightseagreen","red","violetred4","mediumorchid1","blue","royalblue","steelblue1","orange","gold","black")', file=ROUT)
                
                print('y_up_bound = max(c(%s))' % (','.join(["log10(" + str(i) + "+1)" for i in all_genotypes])), file=ROUT)
                print('y_low_bound = min(c(%s))' % (','.join(["log10(" + str(i) + "+1)" for i in all_genotypes])), file=ROUT)
                
                print("pdf(\"%s\")" % (outfile + '.mismatch_profile.pdf'), file=ROUT)
                count=1
                for gt in all_genotypes:
                        if count == 1:
                                print("plot(log10(%s+1),type=\"l\",col=color_code[%d],ylim=c(y_low_bound,y_up_bound),ylab=\"log10(# of mismatch)\",xlab=\"Read position (5\'->3\')\")" % (gt,count), file=ROUT)
                        else:
                                print("lines(log10(%s+1), col=color_code[%d])" % (gt,count), file=ROUT)
                        count += 1
                print("legend(13,y_up_bound,legend=c(%s), fill=color_code, border=color_code, ncol=4)" % (','.join(['"' + i + '"' for i in all_genotypes])), file=ROUT)
                print("dev.off()", file=ROUT)
        def deletionProfile(self,read_length,read_num, outfile, q_cut=30):
                '''
                Calculate deletion profile. 
                Deletion: Deletion from the read (relative to the reference), CIGAR operator 'D'
                '''
                
                DOUT = open(outfile + '.deletion_profile.txt','w')
                ROUT = open(outfile + '.deletion_profile.r','w')
                
                #reading input SAM file
                if self.bam_format:print("Process BAM file ... ", end=' ', file=sys.stderr)
                else:print("Process SAM file ... ", end=' ', file=sys.stderr)
                
                count = 0
                del_postns = collections.defaultdict(int)       #key: position of read. value: deletion times
                #del_sizes = collections.defaultdict(int)       #key: deletion size. value: deletion frequency of this size
                try:
                        while(1):
                                if count >= read_num:
                                        print("Total reads used: " + str(count), file=sys.stderr)
                                        break
                                aligned_read = next(self.samfile)
                                if aligned_read.is_qcfail:continue                      #skip low quanlity                                      
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_unmapped:continue            #skip unmap read
                                if aligned_read.mapq < q_cut: continue
                                if aligned_read.is_reverse:
                                        strand = '-'
                                else:
                                        strand = '+'
                                
                                # skip if read doesn't have deletion
                                read_cigar = aligned_read.cigar
                                if 2 not in [i[0] for i in read_cigar]:         #read contains no deletion
                                        continue

                                # skip partially mapped read
                                read_seq = aligned_read.seq
                                if len(read_seq) != read_length:
                                        continue
                                
                                matched_portion_size = 0
                                for op,value in aligned_read.cigar:
                                        if op == 0: matched_portion_size += value       #match
                                        if op == 4: matched_portion_size += value       #soft clp
                                        if op == 1: matched_portion_size += value   #insertion to read
                                if matched_portion_size != read_length:
                                        continue
                                
                                count += 1
                                del_positions = bam_cigar.fetch_deletion_range(read_cigar)      #[(position, size),(position, size),...]
                                for (p,s) in del_positions:
                                        if strand == '-': p = read_length - p
                                        del_postns[p] += 1
                                        #del_sizes[s] += 1
                                #print aligned_read
                                #print del_positions
                                                
                except StopIteration:
                        print("Total reads used: " + str(count), file=sys.stderr)
                print('\n')
                
                del_count=[]
                print("read_position\tdeletion_count", file=DOUT)
                for k in range(0,read_length):
                        if k in del_postns:
                                print(str(k) + '\t' + str(del_postns[k]), file=DOUT)
                                del_count.append(str(del_postns[k]))
                        else:
                                print(str(k) + '\t0', file=DOUT)
                                del_count.append('0')
                DOUT.close()
                
                print("pdf(\"%s\")" % (outfile + '.deletion_profile.pdf'), file=ROUT)
                print("pos=c(%s)" % ','.join([str(i) for i in range(0,read_length)]), file=ROUT)
                print("value=c(%s)" % ','.join([i for i in del_count]), file=ROUT)              
                print("plot(pos,value,type='b', col='blue',xlab=\"Read position (5\'->3\')\", ylab='Deletion count')", file=ROUT)                                       
                print("dev.off()", file=ROUT)

def print_bits_as_bed( bits ):
        end = 0
        while 1:
                start = bits.next_set( end )
                if start == bits.size: break
                end = bits.next_clear( start )
                print("%d\t%d" % ( start, end ))





























