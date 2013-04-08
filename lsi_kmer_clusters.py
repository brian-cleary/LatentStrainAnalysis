import sys, getopt
import glob, os
from gensim import corpora
import numpy as np
from eigenhashes import Eigenhashes

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
	hashobject = Eigenhashes(inputdir,outputdir)
	corpus_tfidf = corpora.MmCorpus(hashobject.output_path+'initial_corpus.mm')
	lsi = hashobject.train_lsi(corpus_tfidf,num_dims=10,distributed=True)
	lsi.save(hashobject.output_path+'lsi_model.gensim')
	ER = np.transpose(lsi.projection.u*lsi.projection.s)
	C = hashobject.kmer_clusters(ER)
	for c in range(len(C)):
		np.save(hashobject.output_path+str(c)+'.cluster.npy',c)
	np.save(hashobject.output_path+'kmer_cluster_sizes.npy',[len(c) for c in C])