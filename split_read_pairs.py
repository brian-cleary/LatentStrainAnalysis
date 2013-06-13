#!/usr/bin/env python

import sys, getopt
import glob, os
from fastq_reader import Fastq_Reader

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
	FP = glob.glob(os.path.join(inputdir,'*.fastq'))
	FP = [fp for fp in FP if (('.mate1.fastq' not in fp) and ('.mate2.fastq' not in fp) and ('.singleton.fastq' not in fp))]
	FP.sort()
	file_prefix = FP[fr]
	file_prefix = file_prefix[file_prefix.rfind('/')+1:file_prefix.index('.fastq')]
	read_count = hashobject.sort_read_pairs(file_prefix)
	if read_count > 0:
		print file_prefix,'READ COUNT:',str(read_count)
	else:
		print file_prefix,'NO READS'