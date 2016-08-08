#!/usr/bin/env python
'''---------------------------------------------------------------------------------------
Calculate Hexamer frequency
------------------------------------------------------------------------------------------'''

import os,sys
import string
from optparse import OptionParser
import warnings
import string
from qcmodule import FrameKmer

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def file_exist(file):
	try:
   		with open(file) as f: return True
	except IOError as e:
  	 return False


def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",dest="input_read",help="Read sequence in fasta or fastq format. Multiple fasta/fastq files should be separated by ','. For example: read.fq,read2.fa,read3,fa ")
	parser.add_option("-r","--refgenome",action="store",type="string",dest="ref_genome",help="Reference genome sequence in fasta format. Optional")
	parser.add_option("-g","--refgene",action="store",type="string",dest="ref_gene",help="Reference mRNA sequence in fasta format. Optional")
	(options,args)=parser.parse_args()
	
	if not options.input_read:
		parser.print_help()
		sys.exit(0)

	read_table={}
	read_file_names=[]	#base name
	read_file_sum = {}	#sum of hexamer
	
	for read_file in options.input_read.split(','):
		if not file_exist(read_file):
			print >>sys.stderr, read_file, ' does NOT exist!'
			continue	
		print >>sys.stderr, "Calculate hexamer of " + read_file + ' file ...',
		read_table[os.path.basename(read_file)] =  FrameKmer.kmer_freq_file(fastafile = read_file, word_size = 6, step_size = 1, frame = 0)
		read_file_names.append(os.path.basename(read_file))
		read_file_sum[os.path.basename(read_file)] = float(sum(read_table[os.path.basename(read_file)].values()))
		print >>sys.stderr, "Done"	
		
	if options.ref_genome and file_exist(options.ref_genome):
		print >>sys.stderr, "Calculate hexamer of " + options.ref_genome + ' file ...',
		read_table[os.path.basename(options.ref_genome)] = FrameKmer.kmer_freq_file(fastafile = options.ref_genome, word_size = 6, step_size = 1, frame = 0)
		read_file_names.append(os.path.basename(options.ref_genome))
		read_file_sum[os.path.basename(options.ref_genome)] = float(sum(read_table[os.path.basename(options.ref_genome)].values()))
		print >>sys.stderr, "Done."
		
	if options.ref_gene and file_exist(options.ref_gene):
		print >>sys.stderr, "Calculate hexamer of " + options.ref_gene + ' file ...',
		read_table[os.path.basename(options.ref_gene)]= FrameKmer.kmer_freq_file(fastafile = options.ref_gene, word_size = 6, step_size = 1, frame = 0)	
		read_file_names.append(os.path.basename(options.ref_gene))
		read_file_sum[os.path.basename(options.ref_gene)] = float(sum(read_table[os.path.basename(options.ref_gene)].values()))
		print >>sys.stderr, "Done."
	
	print '\n\nHexamer' + '\t' + '\t'.join(read_file_names)
		
	for kmer in FrameKmer.all_possible_kmer(6):
		if 'N' in kmer:continue
		print kmer + '\t',
		try:
			print '\t'.join([str(read_table[name][kmer] / (read_file_sum[name])) for name in read_file_names])
		except:
			print '\t'.join([str(read_table[name][kmer] / (read_file_sum[name]+1)) for name in read_file_names])
		
if __name__ == '__main__':
	main()
