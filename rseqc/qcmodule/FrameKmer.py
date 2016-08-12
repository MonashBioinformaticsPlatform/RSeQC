#!/usr/bin/env python
'''deal with Kmer. DNA sequence should only A, C, G, T. python2.7 or newer'''

#import built-in modules
import os,sys
import numpy
import math
from collections import Counter
import re
import itertools

def word_generator(seq,word_size,step_size,frame=0):
	'''generate DNA word from sequence using word_size and step_size. Frame is 0, 1 or2'''
	for i in xrange(frame,len(seq),step_size):
		word =  seq[i:i+word_size]
		if len(word) == word_size:
			yield word

def seq_generator(fastafile):
	'''DNA sequence only contains A,C,G,T,N. sequence with other characters will be removed'''
	tmpseq=''
	name=''
	DNA_pat = re.compile(r'^[ACGTN]+$')
	for line in open(fastafile,'r'):
		line=line.strip().upper()
		if line.startswith(('#',' ','\n')):continue
		if line.startswith(('>','@')):
			if tmpseq:
				yield [name,tmpseq]
				tmpseq=''
			name = line.split()[0][1:]
		elif DNA_pat.match(line):
			tmpseq += line
	yield [name,tmpseq]
	
def all_possible_kmer(l):
	'''return all possible combinations of A,C,G,T,N. only support A,C,G,T,N. l is length of kmer'''
	for i in itertools.product(['A','C','G','T','N'],repeat=l):
		yield ''.join(i)

def kmer_freq_file (fastafile,word_size,step_size=1,frame=0,min_count=0):
	'''Calculate kmer frequency from fasta file'''
	seq_num = 0
	ret_dict={}
	for n,s in seq_generator(fastafile):
		seq_num += 1
		if seq_num == 1:
			count_table = Counter(word_generator(s,word_size=word_size,step_size=step_size,frame=frame))
		else:
			count_table.update( word_generator(s,word_size=word_size,step_size=step_size,frame=frame) )
	
	#return count_table
	for kmer in all_possible_kmer(word_size):
		if not count_table.has_key(kmer): count_table[kmer]=0
		if count_table[kmer] >= min_count:
			#print kmer + '\t' + str(count_table[kmer])
			if 'N' in kmer:continue
			ret_dict[kmer] = count_table[kmer]
	return ret_dict		

def kmer_freq_seq (seq,word_size,step_size=1,frame=0,min_count=0):
	'''Calculate kmer frequency from DNA sequence. coding. genome is hexamer table calculated
	from coding region and whole genome (as background control)
	'''
	count_table = Counter(word_generator(seq,word_size=word_size,step_size=step_size,frame=frame))
	for kmer in all_possible_kmer(word_size):
		if not count_table.has_key(kmer): count_table[kmer]=0
		if count_table[kmer] >= min_count:
			print kmer + '\t' + str(count_table[kmer])

def kmer_ratio(seq,word_size,step_size,coding,noncoding):
	if len(seq) < word_size:
		return 0
		
	sum_of_log_ratio_0 = 0.0
	sum_of_log_ratio_1 = 0.0
	sum_of_log_ratio_2 = 0.0	
	frame0_count=0.0
	frame1_count=0.0
	frame2_count=0.0
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=0):	
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_0  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_0 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_0 -= 1
		else:
			continue
		frame0_count += 1
	'''	
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=1):
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_1  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_1 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_1 -= 1
		else:
			continue
		frame1_count += 1
	
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=2):
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_2  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_2 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_2 -= 1
		else:
			continue
		frame2_count += 1
	return max(sum_of_log_ratio_0/frame0_count, sum_of_log_ratio_1/frame1_count,sum_of_log_ratio_2/frame2_count)	
	'''
	return sum_of_log_ratio_0/frame0_count
	

