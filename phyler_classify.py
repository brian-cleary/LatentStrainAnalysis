#!/usr/bin/env python

import sys, getopt
import glob,os

# sample the first 10**7 reads
def get_fasta(fp,fo):
	f = open(fp)
	g = open(fo,'w')
	lastlinechar = ''
	writenext = False
	read_count = 0
	for line in f:
		if (line[0] == '@') and (lastlinechar != '+'):
			g.write('>'+line[1:])
			writenext = True
			read_count += 1
		elif writenext:
			g.write(line)
			writenext = False
		lastlinechar = line[0]
		if read_count >= 10**7:
			break
	f.close()
	g.close()
	return read_count

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
	os.system('mkdir '+outputdir+fr)
	FP = glob.glob(os.path.join(inputdir+fr,'*.fastq'))
	read_count = 0
	for fp in FP:
		fileprefix = fp[fp.rfind('/')+1:fp.index('.fastq')]
		fasta_file = outputdir + fr + fileprefix + '.fasta'
		read_count += get_fasta(fp,fasta_file)
	os.system('cat %s*.fasta > %sall.fa' % (outputdir+fr,outputdir+fr))
	os.system('rm '+outputdir+fr+'*.fasta')
	os.system('touch '+outputdir + fr + 'all.count.' + str(read_count))
	os.system('blastall -p blastn -W15 -a1 -e0.01 -m8 -b1 -i %s -d /seq/msctmp/bcleary/src/MetaPhylerV1.25/markers/markers.dna > %s' % (outputdir+fr+'all.fa',outputdir+fr+'all.phyler.blastn'))
	os.system('rm '+outputdir+fr+'all.fa')
	os.system('/seq/msctmp/bcleary/src/MetaPhylerV1.25/metaphylerClassify /seq/msctmp/bcleary/src/MetaPhylerV1.25/markers/markers.blastn.classifier /seq/msctmp/bcleary/src/MetaPhylerV1.25/markers/markers.taxonomy %s > %s' % (outputdir+fr+'all.phyler.blastn',outputdir+fr+'all.phyler.blastn.classification'))
	os.system('rm '+outputdir+fr+'all.phyler.blastn')