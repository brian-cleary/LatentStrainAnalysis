#!/usr/bin/env python

import sys, getopt
import glob, os
from collections import defaultdict

job_name_key = ('#BSUB -J',2)
log_str_keys = {}
log_float_keys = {'Total time': 2,'CPU time': 3,'Max Memory': 3}
# for multiple runs in one log file, this will record only the last
def parse_log_file(fp):
	f = open(fp)
	L = f.readlines()
	f.close()
	R = {}
	try:
		for l in L:
			if l[:len(job_name_key[0])] == job_name_key[0]:
				job = l.strip().split()[job_name_key[1]]
			for k in log_str_keys.keys():
				if l.strip()[:len(k)] == k:
					R[k] = l.strip().split()[log_str_keys[k]]
					break
			for k in log_float_keys.keys():
				if l.strip()[:len(k)] == k:
					R[k] = float(l.strip().split()[log_float_keys[k]])
					break
		return (job,R)
	except Exception,err:
		print str(err)
		return None

help_message = 'usage example: python parse_logs.py -i /project/Logs/ -o /project/log_summary.txt'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:',["inputdir=","outputfile="])
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
		elif opt in ('-o','--outputfile'):
			outfile = arg
	FP = glob.glob(os.path.join(inputdir,'*.out'))
	L = {}
	for fp in FP:
		results = parse_log_file(fp)
		if results != None:
			if results[0] not in L:
				L[results[0]] = defaultdict(list)
			for k,v in results[1].items():
				L[results[0]][k].append(v)
	f = open(outfile,'w')
	for k,v in L.items():
		f.write(k+'\n')
		for kk,vv in v.items():
			vv.sort()
			f.write('\t'+kk+' min: '+str(vv[0])+'\n')
			f.write('\t'+kk+' max: '+str(vv[-1])+'\n')
			f.write('\t'+kk+' median: '+str(vv[len(vv)/2])+'\n')
		f.write('\n')
	f.close()