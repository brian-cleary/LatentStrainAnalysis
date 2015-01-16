#!/usr/bin/env python

import sys, getopt
import glob, os
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
	FP = glob.glob(os.path.join(inputdir,'*.hashq.*'))
	FP = [fp[fp.rfind('/')+1:] for fp in FP]
	FP = list(set([fp[:fp.index('.')] for fp in FP]))
	file_prefix = FP[fr%len(FP)]
	# SUPER DUMB to hardcode the fraction size
	file_fraction = fr/len(FP)
	hashobject = Fastq_Reader(inputdir,outputdir)
	H = hashobject.hash_counts_from_hashq(file_prefix,multi_files_fraction=file_fraction)