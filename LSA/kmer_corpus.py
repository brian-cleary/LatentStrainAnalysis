#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np
from streaming_eigenhashes import StreamingEigenhashes

help_message = 'usage example: python kmer_corpus.py -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:r:',["inputdir=","outputdir=","filerank="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-r','--filerank'):
			fr = int(arg)-1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	hashobject = StreamingEigenhashes(inputdir,outputdir,get_pool=False)
	Kmer_Hash_Count_Files = glob.glob(os.path.join(hashobject.input_path,'*.count.hash'))
	#M = np.load(hashobject.input_path+'column_mask.npy')
	M = []
	hashobject.kmer_corpus_to_disk(Kmer_Hash_Count_Files[fr],mask=M)