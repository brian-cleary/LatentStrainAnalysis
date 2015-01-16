#!/usr/bin/env python

import glob,os
import sys, getopt
import gzip
from fastq_reader import Fastq_Reader

# PAIRED READ FILES ARE ASSUMED TO BE SORTED
def kmer_bins(b,A,pfx,outfile,type):
	if type == 1:
		# use this for readid 1, readid 2 pairs
		def get_id(a):
			return a[:a.index(' ')+2]
	elif type == 2:
		# use this for readid/1, readid/2 pairs
		def get_id(a):
			return a.split()[0]
	else:
		# no known read type treated as singleton
		def get_id(a):
			return a.split()[0]+'*'
	current_id = None
	pair = []
	bins = []
	reads_hashed = 0
	for a in range(len(A)):
		read_id = get_id(A[a])
		if read_id != current_id:
			if (len(bins) > 0) and (read_id[:-1] != current_id[:-1]):
				for rp in pair:
					outfile.write(rp)
					outfile.write(pfx+','.join([str(x) for x in bins]) + ']\n')
					reads_hashed += 1
				pair = []
				bins = []
			current_id = read_id
			pair.append(A[a])
		bins.append(b[a])
	return reads_hashed

help_message = 'usage example: python hash_fastq_reads.py -r 1 -i /project/home/original_reads/ -o /project/home/hashed_reads/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:z',["filerank=","inputdir=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	do_reverse_compliment = True
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
		elif opt in ('-z','--reversecomp'):
			do_reverse_compliment = False
	FP = glob.glob(os.path.join(inputdir,'*.fastq.*'))
	if len(FP) == 0:
		# single file per-sample
		FP = glob.glob(os.path.join(inputdir,'*.fastq'))
	file_prefix = FP[fr]
	file_split = file_prefix[file_prefix.index('.fastq')+6:]
	file_prefix = file_prefix[file_prefix.rfind('/')+1:file_prefix.index('.fastq')]
	hashobject = Fastq_Reader(inputdir,outputdir)
	f = open(hashobject.input_path+file_prefix+'.fastq'+file_split,'r')
	read_type = hashobject.id_type(f)
	g = gzip.open(hashobject.output_path+file_prefix+'.hashq'+file_split+'.gz','wb')
	hashobject.hpfx = hashobject.hpfx + str(hashobject.kmer_size)+','
	A = []
	reads_hashed = 0
	while A != None:
		try:
			A,B = hashobject.generator_to_bins(hashobject.read_generator(f,max_reads=25000,verbose_ids=True),rc=do_reverse_compliment)
			for b in range(len(B)):
				reads_hashed += kmer_bins(B[b],A,hashobject.hpfx,g,read_type)
		except Exception,err:
			pass
			#print str(err)
	f.close()
	g.close()
	print 'total reads hashed:',reads_hashed