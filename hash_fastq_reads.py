#!/usr/bin/env python

import sys, getopt
from fastq_reader import Fastq_Reader

# PAIRED READ FILES ARE ASSUMED TO BE SORTED
def kmer_bins(b,A,pfx,outfile):
	current_id = None
	pair = []
	bins = []
	for a in range(len(A)):
		read_id = A[a].split()[0]
		if read_id != current_id:
			if (len(bins) > 0) and (read_id[:-1] != current_id[:-1]):
				for rp in pair:
					outfile.write(rp)
					outfile.write(pfx+','.join([str(x) for x in bins]) + ']\n')
				pair = []
				bins = []
			pair.append(A[a])
			current_id = read_id
		bins.append(b[a])

help_message = 'usage example: python hash_fastq_reads.py -f SRS013705 -i /project/home/original_reads/ -o /project/home/hashed_reads/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hf:r:i:o:',["fileprefix=","filerank=","inputdir=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-f','--fileprefix'):
			file_prefix = arg
		elif opt in ('-r','--filerank'):
			file_prefix = int(arg) - 1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	if isinstance(file_prefix,int):
		import glob, os
		FP = glob.glob(os.path.join(inputdir,'*.fastq.*'))
		FP = [fp for fp in FP if 'random_kmers' not in fp]
		fp = FP[file_prefix]
		file_split = fp[-3:]
		file_prefix = fp[:-3][fp.rfind('/')+1:]
		file_prefix = file_prefix[:file_prefix.index('.')]
	hashobject = Fastq_Reader(inputdir,outputdir)
	f = open(hashobject.input_path+file_prefix+'.fastq'+file_split,'r')
	g = open(hashobject.output_path+file_prefix+'.hashq'+file_split,'w')
	hashobject.hpfx = hashobject.hpfx + str(hashobject.kmer_size)+','
	A = []
	while A != None:
		try:
			A,B = hashobject.generator_to_bins(hashobject.read_generator(f,max_reads=10000,verbose_ids=True),rc=True)
			for b in range(len(B)):
				kmer_bins(B[b],A,hashobject.hpfx,g)
		except Exception,err:
			print str(err)
	f.close()
	g.close()