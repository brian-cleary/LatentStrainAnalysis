#!/usr/bin/env python

import sys, getopt
import glob, os
from streaming_eigenhashes import StreamingEigenhashes

help_message = 'usage example: python kmer_clusters.py -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/ -p 16'
if __name__ == "__main__":
	singleInstance = False
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:p:s',["inputdir=","outputdir=","numproc"])
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
		elif opt in ('-p','--numproc'):
			num_proc = int(arg)
		elif opt in ('-s','--single'):
			singleInstance = True
	### use -p option for multiprocessing
	num_proc = -1
	###
	hashobject = StreamingEigenhashes(inputdir,outputdir,get_pool=num_proc)
	Kmer_Hash_Count_Files = glob.glob(os.path.join(hashobject.input_path,'*.count.hash.conditioned'))
	hashobject.path_dict = {}
	for i in range(len(Kmer_Hash_Count_Files)):
		hashobject.path_dict[i] = Kmer_Hash_Count_Files[i]
	corpus = hashobject.kmer_corpus_from_disk()
	# This is a hack. Should do a better job chosing num_dims
	lsi = hashobject.train_kmer_lsi(corpus,num_dims=len(hashobject.path_dict)*4/5,single=singleInstance)
	lsi.save(hashobject.output_path+'kmer_lsi.gensim')
	