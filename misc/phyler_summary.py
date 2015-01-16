#!/usr/bin/env python

import sys,getopt
import glob,os
from collections import defaultdict
from operator import itemgetter

def tax_level_summary(FP,inverted=True):
	Classifications = defaultdict(float)
	Samples = {}
	Totals = defaultdict(float)
	for fp in FP:
		L = open(fp).readlines()
		n = fp[:fp.rfind('/')]
		n = n[n.rfind('/')+1:]
		Samples[n] = defaultdict(float)
		for l in L[1:]:
			ls = l.strip().split('\t')
			Classifications[ls[0]] += float(ls[2])
			Samples[n][ls[0]] = ls[1]
			Totals[n] += float(ls[2])
	Classifications = sorted(Classifications.items(),key=itemgetter(1),reverse=True)
	Totals = sorted(Totals.items(),key=itemgetter(1),reverse=True)
	if inverted:
		outlines = ['\t'.join(['Name'] + [k for k,v in Classifications])+'\n']
		outlines.append('\t'.join(['Total Abundance'] + [str(v) for k,v in Classifications])+'\n')
		outlines.append('ID\n')
		for k,v in Totals:
			outlines.append('\t'.join([k] + [str(Samples[k][c[0]]) for c in Classifications])+'\n')
	else:
		outlines = ['\t'.join(['ID']+[t[0] for t in Totals])+'\n']
		outlines.append('\t'.join(['Marker Abundance']+[str(t[1]) for t in Totals])+'\n')
		outlines.append('Name\n')
		for k,v in Classifications:
			outlines.append('\t'.join([k]+[str(Samples[t[0]][k]) for t in Totals])+'\n')
	return outlines

help_message = 'usage example: python phyler_summary.py -i /project/home/phyler/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:',["inputdir="])
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
	FP = glob.glob(inputdir+'*/all.blastn.genus.taxprof.norm')
	f = open(inputdir+'sub-sample_genus_abundance.txt','w')
	f.writelines(tax_level_summary(FP))
	f.close()
	FP = glob.glob(inputdir+'*/all.blastn.family.taxprof.norm')
	f = open(inputdir+'sub-sample_family_abundance.txt','w')
	f.writelines(tax_level_summary(FP))
	f.close()
	FP = glob.glob(inputdir+'*/all.blastn.order.taxprof.norm')
	f = open(inputdir+'sub-sample_order_abundance.txt','w')
	f.writelines(tax_level_summary(FP))
	f.close()