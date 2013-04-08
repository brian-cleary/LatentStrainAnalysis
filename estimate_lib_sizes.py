#!/usr/bin/env python

import sys, getopt
import glob,os
import csv
from fastq_reader import Fastq_Reader

def parse_velvet_log(path):
	f = open(path+'LibLengths')
	l = f.readline()
	f.close()
	try:
		mean = int(l.split()[6][:-1])
		sd = int(l.split()[-1])
		return (mean-sd,mean+sd)
	except:
		return None

help_message = 'usage example: python estimate_lib_sizes.py -i /project/home/original_reads/ -o /project/home/lib_estimates/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:',["inputdir=","outputdir="])
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
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	hashobject = Fastq_Reader(outputdir,outputdir)
	FP = glob.glob(os.path.join(inputdir,'*.fastq.*'))
	FP = set([fp[fp.rfind('/')+1:fp.index('.fastq')] for fp in FP])
	LibSizes = []
	for sample in FP:
		os.system('head -n600000 '+inputdir+sample+'.fastq.aa > '+outputdir+sample+'.fastq')
		sample_reads = hashobject.sort_read_pairs(sample)
		if sample_reads > 0:
			velvetdir = outputdir+sample+'_velvet/'
			os.system('mkdir '+velvetdir)
			os.system('/import/analysis/comp_bio/metagenomics/src/velvet/velveth '+velvetdir+' 31 -fastq -short '+outputdir+sample+'.singleton.fastq -shortPaired '+outputdir+sample+'.pairs.fastq')
			os.system('/import/analysis/comp_bio/metagenomics/src/velvet/velvetg '+velvetdir+' -exp_cov auto | grep -ir "Paired-end library 1 has length:" > '+velvetdir+'LibLengths')
			libsize = parse_velvet_log(velvetdir)
			if libsize != None:
				LibSizes.append((libsize,sample))
		os.system('rm -r '+outputdir+sample+'*')
	LibSizes.sort()
	LibMembers = []
	current_window = (-1,-1)
	current_members = []
	# this is fairly dumb
	for l in LibSizes:
		lib,sample = l
		if lib[0] > current_window[1]:
			if len(current_members) > 0:
				LibMembers.append(current_members)
			current_members = []
			current_window = lib
		current_members.append(sample)
	if len(current_members):
		LibMembers.append(current_members)
	f = open(outputdir+'samples_grouped_by_lib.csv','w')
	writer = csv.writer(f)
	for l in LibMembers:
		writer.writerow(l)
	f.close