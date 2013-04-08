#!/usr/bin/env python

import sys, getopt
import glob, os

help_message = 'usage example: python split_fastq_files.py -i /project/home/original_reads'
if __name__ == "__main__":
	s = 500
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:s:',["inputdir=","splitsize="])
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
		elif opt in ('-s','--splitsize'):
			s = int(arg)
	FP = glob.glob(os.path.join(inputdir,'*.fastq'))
	for fp in FP:
		os.system('split --bytes='+str(s)+'m '+fp+' '+fp+'.')
		os.system('rm '+fp)
		os.system('touch '+fp)