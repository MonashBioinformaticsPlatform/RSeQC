import os,sys
from os.path import abspath,join, getsize
"""
Get bam files from input, input could be:
1) directory that containing one or more bam files
2) plain text file containing one or more bam file paths
3) a single bam file
4) ',' separated multiple bam files

in all cases, the index .bai file(s) should be exist in the same location. eg, if test.bam
exists, test.bam.bai must also exist.
"""

def isbamfile(infile):
	'''check if it is bam file, if it is empty and if the .bai file exists'''
	if os.path.isfile(infile) and infile[-4:].lower() == '.bam':
		if getsize(infile) != 0:
			if os.path.isfile(infile + '.bai'):
				return True
			else:
				print >>sys.stderr, "Warning: %s.bai does not exist! Skip it." % (infile)
				return False
		else:
			print >>sys.stderr, "The size of %s is 0! Skip it." % (infile)
			return False
	else:
		return False	

def get_bam_files (input,printit=False):
	bam_files = []	
	
	#dir
	if os.path.isdir(input):
		for root, directories, files in os.walk(input,followlinks=True):
			full_names = [join(abspath(root), name) for name in files]
			for fn in full_names:
				if isbamfile(fn):
					bam_files.append(fn)
	#single bam file
	elif isbamfile(input):
		bam_files.append(input)
	#plain text file
	elif os.path.isfile(input):
		try:
			for line in open(input):
				line = line.strip()
				if line.startswith('#'):continue
				if isbamfile(line):
					bam_files.append(line)
		except:
			pass
	else:
		tmp = input.split(',')
		if len(tmp) <2: pass
		for i in tmp:
			if isbamfile(i):
				bam_files.append(i)
	
	if printit:
		for i in bam_files:
			print i		
	return bam_files


if __name__ == '__main__':
	get_bam_files(sys.argv[1])