#!/usr/bin/env python

import sys, getopt
import glob, os
import cPickle
import csv
from operator import itemgetter
from collections import defaultdict
from fastq_reader import Fastq_Reader

bact_names_path = '/import/scratch/user/blclea/data/bacterial_seqid_names.cpickle'
vir_names_path = '/import/scratch/user/blclea/data/viral_seqid_names.cpickle'

help_message = 'usage example: python assembly_summary.py -i /project/home/'
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
	hashobject = Fastq_Reader(inputdir+'read_partitions/',inputdir+'read_partitions/')
	f = open(bact_names_path)
	BNames = cPickle.load(f)
	f.close()
	f = open(vir_names_path)
	VNames = cPickle.load(f)
	f.close()
	f = open(inputdir+'lib_estimates/samples_grouped_by_lib.csv')
	reader = csv.reader(f)
	Sample_ids = []
	for row in reader:
		Sample_ids += row
	f.close()
	f_main = open(inputdir+'assembly_alignment_summary.csv','w')
	writer_main = csv.writer(f_main)
	writer_main.writerow(['partition','N50','largest contig','total bp','scaffolds','top bacterial alignment','alignment length','top viral alignment','alignment length'])
	f_bact = open(inputdir+'bacterial_alignments.csv','w')
	writer_bact = csv.writer(f_bact)
	writer_bact.writerow(['partition'] + ['alignment','alignment length']*10)
	f_vir = open(inputdir+'viral_alignments.csv','w')
	writer_vir = csv.writer(f_vir)
	writer_vir.writerow(['partition'] + ['alignment','alignment length']*10)
	f_counts = open(inputdir+'partition_read_counts_by_sample.csv','w')
	writer_counts = csv.writer(f_counts)
	writer_counts.writerow(['partition']+Sample_ids)
	f_genus = open(inputdir+'partition_abundance_genus.csv','w')
	writer_genus = csv.writer(f_genus)
	Genus = defaultdict(int)
	PGenus = {}
	f_family = open(inputdir+'partition_abundance_family.csv','w')
	writer_family = csv.writer(f_family)
	Family = defaultdict(int)
	PFamily = {}
	FP = glob.glob(os.path.join(inputdir+'cluster_vectors/','*.cluster.npy'))
	Clusters = [fp[fp.rfind('/')+1:fp.index('.npy')] for fp in FP]
	for c in Clusters:
		hashobject.assembly_loc = hashobject.output_path+c+'_velvet/'
		try:
			asbly = hashobject.assembly_stats(c)
		except:
			asbly = {}
		try:
			algn = hashobject.alignment_stats(c)
			algn = sorted(algn.iteritems(),key=itemgetter(1),reverse=True)
			if len(algn) == 0:
				algn = [(None,0)]
		except:
			algn = [(None,0)]
		try:
			algn_v = hashobject.alignment_stats(c+'.viral')
			algn_v = sorted(algn_v.iteritems(),key=itemgetter(1),reverse=True)
			if len(algn_v) == 0:
				algn_v = [(None,0)]
		except:
			algn_v = [(None,0)]
		sample_counts = []
		for sample in Sample_ids:
			try:
				sample_counts.append(sum(hashobject.read_counts(c+'.'+sample).values()))
			except:
				sample_counts.append(0)
		writer_main.writerow([c[:c.index('.cluster')],asbly.get('N50',0),asbly.get('largest contig',0),asbly.get('total bp',0),asbly.get('num scaffolds',0),BNames.get(algn[0][0],''),algn[0][1],VNames.get(algn_v[0][0],''),algn_v[0][1]])
		brow = [c[:c.index('.cluster')]]
		for a in algn[:10]:
			brow += [BNames.get(a[0],''),a[1]]
		writer_bact.writerow(brow)
		vrow = [c[:c.index('.cluster')]]
		for a in algn_v[:10]:
			vrow += [VNames.get(a[0],''),a[1]]
		writer_vir.writerow(vrow)
		writer_counts.writerow([c[:c.index('.cluster')]] + sample_counts)
		try:
			A = hashobject.phyler_stats(c,level='genus')
			PGenus[c] = {}
			for a in A:
				Genus[a[0]] += a[2]
				PGenus[c][a[0]] = a[1]
		except:
			PGenus[c] = {}
		try:
			A = hashobject.phyler_stats(c,level='family')
			PFamily[c] = {}
			for a in A:
				Family[a[0]] += a[2]
				PFamily[c][a[0]] = a[1]
		except:
			PFamily[c] = {}
	Genus = sorted(Genus.items(),key=itemgetter(1),reverse=True)
	Genus = [x[0] for x in Genus]
	writer_genus.writerow(['partition']+Genus)
	Family = sorted(Family.items(),key=itemgetter(1),reverse=True)
	Family = [x[0] for x in Family]
	writer_family.writerow(['partition']+Family)
	for c in Clusters:
		writer_genus.writerow([c[:c.index('.cluster')]] + [PGenus[c].get(x,0) for x in Genus])
		writer_family.writerow([c[:c.index('.cluster')]] + [PFamily[c].get(x,0) for x in Family])
	f_main.close()
	f_bact.close()
	f_vir.close()
	f_counts.close()
	f_genus.close()
	f_family.close()