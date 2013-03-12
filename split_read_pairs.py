#!/usr/bin/env python

import sys, getopt
import glob, os

# ASSUMES 4 LINES PER READ, INTERLEAVED SORTING
def split_pairs(i,o):
	f = open(i)
	g1 = open(o+'.1.fastq','w')
	g2 = open(o+'.2.fastq','w')
	firstread = True
	line = f.readline()
	while line != '':
		if firstread:
			g1.write(line)
			for j in range(3):
				g1.write(f.readline())
			firstread = False
		else:
			g2.write(line)
			for j in range(3):
				g2.write(f.readline())
			firstread = True
		line = f.readline()
	f.close()
	g1.close()
	g2.close()


help_message = 'usage example: python split_read_pairs.py -i /project/somefile.pairs.fastq -o /project/somefile'
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
		elif opt in ('-o','--outputdir'):
			outputdir = arg
	split_pairs(inputdir,outputdir)