#!/usr/bin/env python

import sys, getopt
import glob, os
from numpy import load
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
	Kmer_Hash_Count_Files = glob.glob(os.path.join(hashobject.input_path,'*.count.hash'))
	M = hashobject.matrix_from_file_paths(Kmer_Hash_Count_Files)
	M,NZ = hashobject.nonzero_abundances(M)
	M = hashobject.conditioned_nonzeros(M,NZ)
	NZ = None
	# seems to throw segfault with hash size 29 and num_dims > 8
	ER = hashobject.eigenkmers(M,num_dims=8)
	M = None
	C = hashobject.kmer_clusters(ER)
	NZ = load(hashobject.output_path+'nonzero_indices.npy')
	hashobject.save_clusters(C,NZ)

