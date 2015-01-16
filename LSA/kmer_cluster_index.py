#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np
from gensim import models
from streaming_eigenhashes import StreamingEigenhashes

help_message = 'usage example: python kmer_cluster_part.py -r 1 -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:t:',["filerank=","inputdir=","outputdir=","thresh="])
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
		elif opt in ('-r','--filerank'):
			fr = int(arg) - 1
		elif opt in ('-t','--thresh'):
			thresh = float(arg)
	hashobject = StreamingEigenhashes(inputdir,outputdir,get_pool=-1)
	Kmer_Hash_Count_Files = glob.glob(os.path.join(hashobject.input_path,'*.count.hash.conditioned'))
	hashobject.path_dict = {}
	for i in range(len(Kmer_Hash_Count_Files)):
		hashobject.path_dict[i] = Kmer_Hash_Count_Files[i]
	lsi = models.LsiModel.load(hashobject.output_path+'kmer_lsi.gensim')
	hashobject.cluster_thresh = thresh
	Index = hashobject.lsi_cluster_index(lsi)
	np.save(hashobject.output_path+'cluster_index.npy',Index)
	print 'cluster index has shape:',Index.shape
	f = open(hashobject.output_path+'numClusters.txt','w')
	f.write('%d\n' % Index.shape[0])
	f.close()