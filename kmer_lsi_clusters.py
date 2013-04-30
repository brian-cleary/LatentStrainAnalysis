#!/usr/bin/env python

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
	lsi = hashobject.train_kmer_lsi(corpus,num_dims=25)
	lsi.save(hashobject.output_path+'kmer_lsi.gensim')
	C = hashobject.lsi_kmer_clusters(lsi)
	for i in range(len(C)):
		np.save(hashobject.output_path+str(i)+'.cluster.npy',C[i])
	np.save(hashobject.output_path+'kmer_cluster_sizes.npy',[len(c) for c in C])
	GW = np.load(hashobject.output_path+'global_weights.npy')
	CP = []
	global_weight_sum = GW.sum(dtype=np.float64)
	for c in C:
		CP.append(GW[c].sum(dtype=np.float64)/global_weight_sum)
	np.save(hashobject.output_path+'cluster_probs.npy',CP)
	I = dict([(k,iter(v)) for k,v in enumerate(C)])
	Ix = defaultdict(list)
	for k,v in I.items():
		Ix[v.next()].append(k)
	f = open(hashobject.output_path+'cluster_cols.txt','w')
	for i,x in enumerate(GW):
		if i in Ix:
			for k in Ix[i]:
				f.write('%d\t%s\t%f\n' % (i,k,x))
				try:
					Ix[I[k].next()].append(k)
				except:
					pass
			del Ix[i]
	f.close()
	