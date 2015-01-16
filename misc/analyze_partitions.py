#!/usr/bin/env python

import sys, getopt
import glob, os
import csv
from fastq_reader import Fastq_Reader

help_message = 'usage example: python kmer_clusters.py -r 2 -i /project/home/read_partitions/ -o /project/home/read_partitions/ -l /project/home/lib_estimates/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:l:',["filerank=","inputdir=","outputdir=","libdir="])
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
		elif opt in ('-l','--libdir'):
			libdir = arg
			if libdir[-1] != '/':
				libdir += '/'
	hashobject = Fastq_Reader(inputdir,outputdir)
	Read_Partitions = glob.glob(os.path.join(hashobject.input_path,'*.fastq'))
	Read_Partitions = [fp for fp in Read_Partitions if ('.pairs.' not in fp) and ('.singleton.' not in fp)]
	Read_Partitions = list(set([fp[fp.rfind('/')+1:fp.index('.cluster')+8] for fp in Read_Partitions]))
	Processed_Partitions = glob.glob(os.path.join(hashobject.output_path,'*.cluster_velvet/contigs.fa'))
	Processed_Partitions = [fp[len(hashobject.output_path):fp.index('.cluster')+8] for fp in Processed_Partitions]
	rp = Read_Partitions[fr]
	if rp not in Processed_Partitions:
		f = open(libdir+'samples_grouped_by_lib.csv')
		reader = csv.reader(f)
		hashobject.sample_library = {}
		i = 0
		for row in reader:
			for sample in row:
				hashobject.sample_library[sample] = i
			i += 1
		f.close()
		hashobject.num_libs = i
		hashobject.process_cluster(rp)