import collections
import sys

start_coden = ['ATG']
stop_coden = ["TAG","TAA","TGA"]
def _reverse_comp(seq):
	swap = {"A":"T", "T":"A", "C":"G", "G":"C","N":"N","X":"X"}
	tmp = "".join(swap[b] for b in seq)
	return tmp[::-1]
  	
def longest_orf(seq,strandness,sc=None,tc=None):
	'''find the longest ORF in input mRNA sequence. strand=0 only search '+' strand, strand=1, only 
	search '-' strand, strand=2, search both strand'''
	orf_ranges=collections.defaultdict(list)
	dna_seq = seq.upper()
	possible_orf={}	#[orf-st, orf_end] ==>size
	largest_orf=""
	
	if sc is not None:
		start_coden = sc.strip(',').split(',')
		if len(start_coden)==0:
			print >>sys.stderr,"Unkown start codon"
			sys.exit(1)
		else:
			for cd in start_coden:
				if len(cd) != 3:
					print >>sys.stderr,"Unkown start codon" + str(cd)
					sys.exit(1)
	if tc is not None:
		stop_coden = tc.strip(',').split(',')
		if len(stop_coden)==0:
			print >>sys.stderr,"Unkown stop codon"
			sys.exit(1)
		else:
			for cd in stop_coden:
				if len(cd) != 3:
					print >>sys.stderr,"Unkown stop codon" + str(cd)
					sys.exit(1)	

	strand = strandness
	start_pos = []
	end_pos=[]
	orf_ranges.clear()
	if strand == '-':
		dna_seq=_reverse_comp(dna_seq)
	for sc in start_coden:
		start_found = dna_seq.find(sc)
		while start_found >-1:
			start_pos.append(start_found)
			start_found = dna_seq.find(sc,start_found +1)
	for sc in stop_coden:
		end_found = dna_seq.find(sc)
		while end_found >-1:
			end_pos.append(end_found)
			end_found = dna_seq.find(sc,end_found+1)
	
	for st in start_pos:
		for end in end_pos:
			if end <= st: continue
			if (end - st) % 3==0:
				orf_ranges[st].append(end)
	for k in sorted(orf_ranges):
		possible_orf[str(k) + '\t' + str(min(orf_ranges[k])) + '\t' + strand] = min(orf_ranges[k]) - k
		
		
	for k,v in possible_orf.items():
		if v == max(possible_orf.values()):
			#print "#" + k
			fields=k.split()
			largest_orf = dna_seq[int(fields[0]):int(fields[1])]
	return largest_orf		#could be None, '' or DNA sequencee
	
				
def longest_orf_bed(seq,bedline,sc=None,tc=None):
	'''find the longest ORF in input mRNA sequence. strand=0 only search '+' strand, strand=1, only 
	search '-' strand, strand=2, search both strand'''
	orf_ranges=collections.defaultdict(list)
	dna_seq = seq.upper()
	possible_orf={}	#[orf-st, orf_end] ==>size
	largest_orf=""

	if sc is not None:
		start_coden = sc.strip(',').split(',')
		if len(start_coden)==0:
			print >>sys.stderr,"Unkown start codon"
			sys.exit(1)
		else:
			for cd in start_coden:
				if len(cd) != 3:
					print >>sys.stderr,"Unkown start codon" + str(cd)
					sys.exit(1)
	if tc is not None:
		stop_coden = tc.strip(',').split(',')
		if len(stop_coden)==0:
			print >>sys.stderr,"Unkown stop codon"
			sys.exit(1)
		else:
			for cd in stop_coden:
				if len(cd) != 3:
					print >>sys.stderr,"Unkown stop codon" + str(cd)
					sys.exit(1)	
	
	fields = bedline.split()
	txStart = int(fields[1])
	exon_sizes = [int(i) for i in fields[10].rstrip(',\n').split(',')]
	exon_starts = map(int,fields[11].rstrip(',').split(','))
	exon_starts = map((lambda x: x + txStart),exon_starts)
	exon_ends = map(int,fields[10].rstrip(',').split(','))
	exon_ends = map((lambda x,y:x+y),exon_starts,exon_ends)
	strand = fields[5]
	
	start_pos = []
	end_pos=[]
	orf_ranges.clear()
	if strand == '-':
		dna_seq=_reverse_comp(dna_seq)
	for sc in start_coden:
		start_found = dna_seq.find(sc)
		while start_found >-1:
			start_pos.append(start_found)
			start_found = dna_seq.find(sc,start_found +1)
	for sc in stop_coden:
		end_found = dna_seq.find(sc)
		while end_found >-1:
			end_pos.append(end_found)
			end_found = dna_seq.find(sc,end_found+1)
	
	for st in start_pos:
		for end in end_pos:
			if end <= st: continue
			if (end - st) % 3==0:
				orf_ranges[st].append(end)
	for k in sorted(orf_ranges):
		possible_orf[str(k) + '\t' + str(min(orf_ranges[k])) + '\t' + strand] = min(orf_ranges[k]) - k
	
	for k,v in possible_orf.items():
		if v == max(possible_orf.values()):
			#print k + '\t' + str(v)
			col=k.split()
			cds_st = int(col[0])
			cds_end = int(col[1])
			strand = col[2]
			
			if strand =='+':
				#determine CDS start site
				for exon_size, exon_st in zip(exon_sizes,exon_starts):
					if cds_st > exon_size: 
						cds_st = cds_st - exon_size
						continue
					else:
						cds_pos1 = exon_st + cds_st
						break
				#determine CDS end site
				for exon_size, exon_st in zip(exon_sizes,exon_starts):
					if cds_end > exon_size:
						cds_end = cds_end - exon_size
						continue
					else:
						cds_pos2 = exon_st + cds_end + 3
						break
				fields[6] = str(cds_pos1)
				fields[7] = str(cds_pos2)
				return '\t'.join(fields)

			if strand =='-':
				exon_sizes = exon_sizes[::-1]
				exon_ends = exon_ends[::-1]
				for exon_size, exon_end in zip(exon_sizes,exon_ends):
					if cds_st > exon_size: 
						cds_st = cds_st - exon_size
						continue
					else:
						cds_pos1 = exon_end - cds_st
						break
				for exon_size, exon_end in zip(exon_sizes,exon_ends):
					if cds_end > exon_size:
						cds_end = cds_end - exon_size
						continue
					else:
						cds_pos2 = exon_end - cds_end - 3
						break
				fields[7] = str(cds_pos1)
				fields[6] = str(cds_pos2)
				return '\t'.join(fields)
