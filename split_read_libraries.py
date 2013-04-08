#!/usr/bin/env python

import sys, getopt

if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:',["inputfile=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-i','--inputfile'):
			inputfile = arg
		elif opt in ('-o','--outputdir'):
			outputdir = arg + '.'
	f = open(inputfile)
	F = {}
	line = f.readline()
	while line != '':
		if line[0] == '@':
			# assuming illumina style
			lib_id = line[1:line.index('.')]
			if lib_id not in F:
				F[lib_id] = open(outputdir+lib_id+'.fastq','w')
			F[lib_id].write(line)
			for _ in range(3):
				F[lib_id].write(f.readline())
		line = f.readline()
	f.close()
	for f in F.values():
		f.close()