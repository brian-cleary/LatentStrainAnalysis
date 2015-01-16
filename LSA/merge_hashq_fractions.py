#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np
from fastq_reader import Fastq_Reader

help_message = 'usage example: python merge_hashq_files.py -r 3 -i /project/home/hashed_reads/ -o /project/home/hashed_reads/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["filerank=","inputdir=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-r','--filerank'):
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
	file_prefix = FP[fr]
	hashobject = Fastq_Reader(inputdir,outputdir)
	H = hashobject.merge_count_fractions(file_prefix)
	H = np.array(H,dtype=np.uint16)
	nz = np.nonzero(H)[0]
	np.save(hashobject.output_path+file_prefix+'.nonzero.npy',nz)
	print 'sample %s has %d nonzero elements and %d total observed kmers' % (file_prefix,len(nz),H.sum())