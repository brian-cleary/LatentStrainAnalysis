#!/usr/bin/env python

import sys, getopt
import glob, os
import subprocess
import numpy as np
from collections import defaultdict
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
	# This is a hack. Should do a better job chosing num_dims
	lsi = hashobject.train_kmer_lsi(corpus,num_dims=len(hashobject.path_dict)*4/5)
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
	Y = np.memmap(hashobject.output_path+'cluster_vals.npy',dtype=np.float32,mode='w+',shape=GW.shape)
	Y[:] = GW
	del Y
	del GW
	X = np.memmap(hashobject.output_path+'cluster_cols.npy',dtype=np.int16,mode='w+',shape=(2**hashobject.hash_size,5))
	I = dict([(k+1,iter(v)) for k,v in enumerate(C)])
	Ix = defaultdict(list)
	for k,v in I.items():
		Ix[v.next()].append(k)
	# should maybe flush at some point
	while len(Ix)>0:
		minx = min(Ix.keys())
		X[minx,:len(Ix[minx])] = Ix[minx]
		for k in Ix[minx]:
			try:
				Ix[I[k].next()].append(k)
			except:
				pass
		del Ix[minx]
	del X
	