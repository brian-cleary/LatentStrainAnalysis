#!/usr/bin/env python

import sys, getopt
import glob, os
from gensim import corpora
import numpy as np
from streaming_eigenhashes import StreamingEigenhashes

help_message = 'usage example: python kmer_clusters.py -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:',["inputdir=","outputdir="])
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
	hashobject = StreamingEigenhashes(inputdir,outputdir,get_pool=False)
	Kmer_Hash_Count_Files = glob.glob(os.path.join(hashobject.input_path,'*.nonzero.npy'))
	corpus_generator = hashobject.corpus_idf_from_hash_paths(Kmer_Hash_Count_Files)
	hashobject.train_tfidf(corpus_generator)
	np.save(hashobject.output_path+'global_weights.npy',hashobject.global_weights)