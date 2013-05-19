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
	FP = glob.glob(os.path.join(inputdir,'*.fastq.*'))
	# run on just a portion of the files
	FP = [fp for fp in FP if fp[-2:] in ('ac','ae')]
	FP = [fp[fp.rfind('/')+1:] for fp in FP]
	FP = list(set([fp[:fp.index('.')] for fp in FP]))
	fileprefix = FP[fr]
	FP = glob.glob(os.path.join(outputdir,fileprefix+'.*.classification'))
	os.system('/import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/combine %s > %s' % (' '.join(FP),outputdir+fileprefix+'.blastn.classification'))
	os.system('rm '+' '.join(FP))
	os.system('/import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/taxprof 0.9 %s %s /import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/markers/tid2name.tab' % (outputdir+fileprefix+'.blastn.classification',outputdir+fileprefix+'.blastn'))