#!/usr/bin/env python

import sys, getopt
import glob,os

help_message = 'usage example: python merge_partition_parts.py -i /project/home/read_partitions'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:',["inputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
	FP = glob.glob(os.path.join(inputdir,'*.fastq.*'))
	FP = set([fp[fp.rfind('/')+1:fp.index('.fastq')] for fp in FP])
	for group in FP:
		gp = glob.glob(os.path.join(inputdir,group+'.fastq.*'))
		os.system('cat '+' '.join(gp)+' > '+inputdir+group+'.fastq')
		os.system('rm '+' '.join(gp))