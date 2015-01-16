#!/usr/bin/env python

import sys, getopt
import glob,os
from fastq_reader import Fastq_Reader

help_message = 'usage example: python merge_read_clusters.py -r 1 -i /project/home/cluster_vectors/intermediate_clusters/ -o /project/home/read_partitions/ -H /project/home/hashed_reads/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:H:',["filerank=","inputdir=","outputdir=","headers"])
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
		elif opt in ('-H',"--headers"):
			headerdir = arg
			if headerdir[-1] != '/':
				headerdir += '/'
	hashobject = Fastq_Reader(inputdir,outputdir)
	FP = glob.glob(os.path.join(hashobject.input_path,'*.*'))
	FP = list(set([fp[fp.rfind('/')+1:fp.rfind('.')] for fp in FP]))
	file_group = FP[fr]
	FP = glob.glob(os.path.join(headerdir,'*.hashq.*'))
	originating_header = FP[int(file_group)]
	originating_header = originating_header[originating_header.rfind('/')+1:originating_header.index('.hashq')]
	hashobject.fastq_from_intermediate_output(file_group,originating_header)