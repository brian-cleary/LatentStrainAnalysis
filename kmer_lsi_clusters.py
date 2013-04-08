import sys, getopt
import glob, os
import subprocess
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
	hashobject = StreamingEigenhashes(inputdir,outputdir)
	Kmer_Hash_Count_Files = glob.glob(os.path.join(hashobject.input_path,'*.count.hash.conditioned'))
	hashobject.path_dict = {}
	for i in range(len(Kmer_Hash_Count_Files)):
		hashobject.path_dict[i] = Kmer_Hash_Count_Files[i]
	corpus = hashobject.kmer_corpus_from_disk()
	lsi = hashobject.train_kmer_lsi(corpus,num_dims=15)
	lsi.save(hashobject.output_path+'kmer_lsi.gensim')
	C = hashobject.lsi_kmer_clusters(lsi)
	for i in range(len(C)):
		np.save(hashobject.output_path+str(i)+'.cluster.npy',C[i])
	np.save(hashobject.output_path+'kmer_cluster_sizes.npy',[len(c) for c in C])