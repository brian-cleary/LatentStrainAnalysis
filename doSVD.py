#!/usr/bin/env python

import sys, getopt
import glob, os

help_message = 'usage example: python doSVD.py -i /project/home/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:',["inputdir="])
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
	os.system('mkdir %scluster_vectors' % (inputdir))
	os.system('python create_jobs.py -j GlobalWeights -i %s' % (inputdir))
	os.system('bsub < %sGlobalWeights_Job.q' % (inputdir))
	os.system('python create_jobs.py -j KmerCorpus -i %s' % (inputdir))
	os.system('bsub -w "done(GlobalWeights)" < %sKmerCorpus_ArrayJob.q' % (inputdir))
	os.system('python create_jobs.py -j LSIKmerClusters -i %s' % (inputdir))
	os.system('bsub -w "done(KmerCorpus)" < %sLSIKmerClusters_Job.q' % (inputdir))
	