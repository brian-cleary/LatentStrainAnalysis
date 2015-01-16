#!/usr/bin/env python

import sys, getopt
import glob, os

help_message = 'usage example: python doHash.py -i /project/home/ -k 49 -s 29'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:k:s:',["inputdir=","kmersize=","hashsize="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-k','--kmersize'):
			k = arg
		elif opt in ('-s','--hashsize'):
			s = arg
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
	os.system('mkdir '+inputdir+'Logs')
	os.system('mkdir '+inputdir+'hashed_reads')
	os.system('python create_jobs.py -j CreateHash -i '+inputdir)
	os.system('bsub < '+inputdir+'CreateHash_Job.q')
	os.system('rm '+inputdir+'original_reads/random_kmers.fastq')
	if len(glob.glob(os.path.join(inputdir+'hashed_reads','Wheels.txt'))) == 1:
		os.system('python create_jobs.py -j HashReads -i '+inputdir)
		os.system('bsub < '+inputdir+'HashReads_ArrayJob.q')
		os.system('python create_jobs.py -j MergeHash -i '+inputdir)
		os.system('bsub -w "done(HashReads)" < '+inputdir+'MergeHash_ArrayJob.q')
		os.system('python create_jobs.py -j CombineFractions -i '+inputdir)
		os.system('bsub -w "done(MergeHash)" < '+inputdir+'CombineFractions_ArrayJob.q')