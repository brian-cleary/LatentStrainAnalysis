#!/usr/bin/env python

import sys, getopt
import glob,os

def get_fasta(fp,fo):
	f = open(fp)
	g = open(fo,'w')
	lastlinechar = ''
	writenext = False
	for line in f:
		if (line[0] == '@') and (lastlinechar != '+'):
			g.write('>'+line[1:])
			writenext = True
		elif writenext:
			g.write(line)
			writenext = False
		lastlinechar = line[0]
	f.close()
	g.close()

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
	FP.sort()
	fileprefix = FP[fr]
	fileprefix = fileprefix[fileprefix.rfind('/')+1:fileprefix.index('.fastq')]+fileprefix[-3:]
	fasta_file = outputdir + fileprefix + '.fasta'
	get_fasta(FP[fr],fasta_file)
	os.system('blastall -p blastn -W15 -a1 -e0.01 -m8 -b1 -i %s -d /import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/markers/markers.dna > %s' % (fasta_file,outputdir+fileprefix+'.phyler.blastn'))
	os.system('rm '+fasta_file)
	os.system('/import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/metaphylerClassify /import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/markers/markers.blastn.classifier /import/analysis/comp_bio/metagenomics/src/MetaPhylerV1.25/markers/markers.taxonomy %s > %s' % (outputdir+fileprefix+'.phyler.blastn',outputdir+fileprefix+'.phyler.blastn.classification'))
	os.system('rm '+outputdir+fileprefix+'.phyler.blastn')