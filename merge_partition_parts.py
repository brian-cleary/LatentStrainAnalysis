#!/usr/bin/env python

import sys, getopt
import glob,os

help_message = 'usage example: python merge_partition_parts.py -r 1 -i /project/home/cluster_vectors/ -o read_partitions/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["inputdir="])
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
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	FP = glob.glob(os.path.join(inputdir,'*.fastq.*'))
	FP = list(set([fp[fp.rfind('/')+1:fp.index('.fastq')] for fp in FP]))
	FP.sort()
	group = FP[fr]
	gp = glob.glob(os.path.join(inputdir,group+'.fastq.*'))
	os.system('cat '+' '.join(gp)+' > '+outputdir+group+'.fastq')
	os.system('touch %s.empty' % (gp[0]))
	os.system('rm '+' '.join(gp))