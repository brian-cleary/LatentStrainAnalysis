#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np
from fastq_reader import Fastq_Reader


help_message = 'usage example: python intermediate_read_clusters.py -r 1 -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["--filerank=","inputdir=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-r',"--filerank"):
			fr = int(arg)-1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	hashobject = Fastq_Reader(inputdir,outputdir)
	Hashq_Files = glob.glob(os.path.join(hashobject.input_path,'*.hashq.*'))
	hashobject.infile = Hashq_Files[fr]
	hashobject.outfile = hashobject.output_path + 'intermediate_clusters/' + str(fr)
	hashobject.global_weights = np.load(hashobject.output_path + 'global_weights.npy')
	global_weight_sum = hashobject.global_weights.sum(dtype=np.float64)
	Cluster_Files = glob.glob(os.path.join(hashobject.output_path,'*.cluster.npy'))
	Cluster_Files = [(int(cf[cf.rfind('/')+1:cf.index('.')]),cf) for cf in Cluster_Files]
	cluster_sizes = np.load(hashobject.output_path+'kmer_cluster_sizes.npy')
	total_set_size = 0
	cluster_weights = []
	cluster_keys = []
	outpart = 0
	for ci,cf in Cluster_Files:
		# ignore super clusters and super small clusters
		if cluster_sizes[ci] < 0.2*2**hashobject.hash_size:
			cw = np.load(cf)
			cw_sum_prob = hashobject.global_weights[cw].sum(dtype=np.float64)/global_weight_sum
			if cw_sum_prob > 0.00002:
				cluster_weights.append((set(cw),cw_sum_prob))
				cluster_keys.append(cf[cf.rfind('/')+1:cf.rfind('.')])
				total_set_size += len(cw)
				if total_set_size > 50*10**6:
					hashobject.membership_generator(cluster_weights,cluster_keys,outpart)
					cluster_weights = []
					cluster_keys = []
					total_set_size = 0
					outpart += 1
	if len(cluster_weights) > 0:
		hashobject.membership_generator(cluster_weights,cluster_keys,outpart)