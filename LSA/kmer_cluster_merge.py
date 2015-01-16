#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np

help_message = 'usage example: python kmer_cluster_merge.py -r 1 -i /project/home/cluster_vectors/ -o /project/home/cluster_vectors/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["filerank=","inputdir=","outputdir="])
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
	FP = glob.glob(os.path.join(inputdir+str(fr),'*.npy'))
	if len(FP) > 0:
		C = np.empty((10**9,),dtype=np.uint64)
		ci = 0
		for fp in FP:
			c = np.load(fp)
			C[ci:ci+c.shape[0]] = c
			ci += c.shape[0]
		np.save(outputdir+str(fr)+'.cluster.npy',C[:ci])
		os.system('rm -r %s%d/' % (outputdir,fr))
	else:
		print 'NO FILES:',fr