#!/usr/bin/env python

import sys, getopt
import glob,os

help_message = 'usage example: python read_phyler.py -r 1 -i /project/home/original_reads/ -o /project/home/phyler/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["inputdir="])
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
	fr = str(fr) + '/'
	### this can be done in phyler_classify
	os.system('/seq/msctmp/bcleary/src/MetaPhylerV1.25/taxprof 0.9 %s %s /seq/msctmp/bcleary/src/MetaPhylerV1.25/markers/tid2name.tab' % (outputdir+fr+'all.phyler.blastn.classification',outputdir+fr+'all.blastn'))
	FP = glob.glob(os.path.join(outputdir+fr,'*.count.*'))
	total = sum([int(fp[fp.rfind('.')+1:]) for fp in FP])
	correction_factor = 100.
	FP = glob.glob(os.path.join(outputdir+fr,'*.taxprof'))
	for fp in FP:
		f = open(fp)
		g = open(fp+'.norm','w')
		g.write(f.readline())
		for line in f:
			ls = line.strip().split('\t')
			g.write('%s\t%s\t%f\n' % (ls[0],ls[1],int(ls[2])*correction_factor/total))
		f.close()
		g.close()