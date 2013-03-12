#!/usr/bin/env python

import sys, getopt
from fastq_reader import Fastq_Reader

help_message = 'usage example: python merge_hashq_files.py -r 3 -i /project/home/hashed_reads/ -o /project/home/hashed_reads/'
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
		elif opt in ('-r','--filerank'):
			fr = int(arg) - 1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	import glob, os
	FP = glob.glob(os.path.join(inputdir,'*.hashq.*'))
	FP = list(set([fp[:-3] for fp in FP]))
	fp = FP[fr]
	file_prefix = fp[fp.rfind('/')+1:]
	file_prefix = file_prefix[:file_prefix.index('.')]
	hashobject = Fastq_Reader(inputdir,outputdir)
	H = hashobject.hash_counts_from_hashq(file_prefix,multi_files=True)