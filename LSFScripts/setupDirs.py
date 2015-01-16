#!/usr/bin/env python

import sys,getopt,os

SplitInput_string = """#!/bin/bash
#BSUB -J SplitInput[1-%numSamples%]
#BSUB -o Logs/SplitInput-Out-%I.out
#BSUB -e Logs/SplitInput-Err-%I.err
#BSUB -q week
#BSUB -W 23:58
echo Date: `date`
t1=`date +%s`
sleep ${LSB_JOBINDEX}
python LSFScripts/array_merge.py -r ${LSB_JOBINDEX} -i %input% -o original_reads/
[ $? -eq 0 ] || echo 'JOB FAILURE: $?'
echo Date: `date`
t2=`date +%s`
tdiff=`echo 'scale=3;('$t2'-'$t1')/3600' | bc`
echo 'Total time:  '$tdiff' hours'
"""

help_message = "usage example: python setupDirs.py -i /path/to/reads/ -n numberOfSamples"
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:n:',["inputdir="])
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
		elif opt in ('-n'):
			n = arg
	for dir in ['Logs','original_reads','hashed_reads','cluster_vectors','read_partitions']:
		os.system('mkdir %s' % (dir))
	f = open('LSFScripts/SplitInput_ArrayJob.q','w')
	f.write(SplitInput_string.replace('%numSamples%',n).replace('%input%',inputdir))
	f.close()
