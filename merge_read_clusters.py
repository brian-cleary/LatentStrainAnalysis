#!/usr/bin/env python

import sys, getopt
import glob,os
from fastq_reader import Fastq_Reader

help_message = 'usage example: python merge_read_clusters.py -i /project/home/cluster_vectors/intermediate_clusters/ -o /project/home/read_partitions/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["--filerank=","inputdir=","outputdir="])
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
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	hashobject = Fastq_Reader(inputdir,outputdir)
	FP = glob.glob(os.path.join(hashobject.input_path,'*.*'))
	FP = set([fp[fp.rfind('/')+1:fp.rfind('.')] for fp in FP])
	for file_group in FP:
		hashobject.fastq_from_intermediate_output(file_group)