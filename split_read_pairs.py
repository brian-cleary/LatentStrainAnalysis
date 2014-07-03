#!/usr/bin/env python

import sys, getopt
import glob, os
from fastq_reader import Fastq_Reader

def get_read_block(f,n=10000):
	R = {}
	line = f.readline()
	last_id = None
	read_lines = []
	while (len(R) < n) and (line != ''):
		while (line != '') and not ((line[0] == '@') and ((line[-3:] == '/1\n') or (line[-3:] == '/2\n'))):
			read_lines.append(line)
			line = f.readline()
		if (len(read_lines) == 3) and (len(read_lines[0]) == len(read_lines[2])):
			R[last_id] = read_lines
		read_lines = []
		last_id = line
		line = f.readline()
	return R

def fix_read_pairs(fp):
	f1 = open(fp+'1.fastq')
	f2 = open(fp+'2.fastq')
	g1 = open(fp+'1.fastq.tmp','w')
	g2 = open(fp+'2.fastq.tmp','w')
	R1 = [None]
	R2 = [None]
	while (len(R1) > 0) and (len(R2) > 0):
		R1 = get_read_block(f1)
		R2 = get_read_block(f2)
		for k,v in R1.iteritems():
			read_id = k.strip()[:-1]
			if read_id+'2\n' in R2:
				g1.write(k + ''.join(v))
				g2.write(read_id+'2\n' + ''.join(R2[read_id+'2\n']))
	g1.close()
	g2.close()

help_message = 'usage example: python split_read_pairs.py -r 1 -i /project/somefile.pairs.fastq -o /project/somefile'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["filerank=","inputdir=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-r',"--filerank"):
			fr = int(arg)-1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
		elif opt in ('-o','--outputdir'):
			outputdir = arg
	hashobject = Fastq_Reader(inputdir,outputdir)
	fr = str(fr) + '/'
	FP = glob.glob(os.path.join(inputdir+fr,'*.fastq'))
	FP = [fp for fp in FP if (('.mate1.fastq' not in fp) and ('.mate2.fastq' not in fp) and ('.singleton.fastq' not in fp))]
	for file_prefix in FP:
		file_prefix = fr + file_prefix[file_prefix.rfind('/')+1:file_prefix.index('.fastq')]
		read_count = hashobject.sort_read_pairs(file_prefix)
		if read_count > 0:
			print file_prefix,'READ COUNT:',str(read_count)
		else:
			print file_prefix,'NO READS'
	FP = glob.glob(os.path.join(inputdir+fr,'*.mate1.fastq'))
	for fp in FP:
		base_fp = fp[:fp.index('1.fastq')]
		fix_read_pairs(base_fp)
		if (os.stat(base_fp+'1.fastq.tmp').st_size > .9*os.stat(base_fp+'1.fastq').st_size) and (os.stat(base_fp+'2.fastq.tmp').st_size > .9*os.stat(base_fp+'2.fastq').st_size):
			os.system('mv %s %s' % (base_fp+'1.fastq.tmp',base_fp+'1.fastq'))
			os.system('mv %s %s' % (base_fp+'2.fastq.tmp',base_fp+'2.fastq'))
		else:
			os.system('rm ' + base_fp+'1.fastq.tmp')
			os.system('rm ' + base_fp+'2.fastq.tmp')
			print 'FAILURE FIXING READS',fp