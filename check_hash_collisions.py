#!/usr/bin/env python

import sys, getopt
import glob,os
from fastq_reader import Fastq_Reader

help_message = 'usage example: python check_hash_collisions.py -r 2 -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/'
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
		elif opt in ('-r',"--filerank"):
			fr = int(arg)
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	hashobject = Fastq_Reader(inputdir,outputdir)
	HashFiles = glob.glob(os.path.join(hashobject.input_path,'*.hashq.*'))
	hashobject.infile = HashFiles[fr]
	t,n,H = hashobject.collision_report()
	print 'total hashed kmers:',t
	print 'total pairwise collisions:',n
	print 'collision histogram:',H