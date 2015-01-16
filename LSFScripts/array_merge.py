#!/usr/bin/env python

import glob,os
import sys, getopt

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
		elif opt in ('-r','--filerank'):
			fr = int(arg) - 1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	FP = glob.glob(os.path.join(inputdir,'*.fastq.1'))
	FP.sort()
	fp = FP[fr]
	p1 = fp
	p2 = fp[:-1] + '2'
	s = fp[:fp.index('.fastq')] + '.single.fastq.1'
	o = outputdir + fp[fp.rfind('/')+1:fp.index('.fastq')]
	os.system('python LSFScripts/merge_and_split_pair_files.py -1 %s -2 %s -s %s -o %s' % (p1,p2,s,o))
	os.system('python LSFScripts/merge_and_split_pair_files.py -s %s -o %s' % (s[:-1] + '2',o))
