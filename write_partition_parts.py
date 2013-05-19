#!/usr/bin/env python

import glob,os
import sys,getopt
import gzip
import numpy as np
from bitarray import bitarray
from collections import defaultdict
from fastq_reader import Fastq_Reader

def get_col(fc,fv,x):
	fc.seek(x*10)
	fv.seek(x*4)
	col = np.fromfile(fc,dtype=np.int16,count=5)
	val = np.fromfile(fv,dtype=np.float32,count=1)
	nz = np.nonzero(col)[0]
	return col[nz]-1,val

def max_log_lik_ratio(s,bkg,h1_prob=0.8,thresh=3.84):
	LLR = [(None,None)]
	read_match_sum = s[-1]
	del s[-1]
	v1 = read_match_sum*h1_prob*(1-h1_prob)
	m1 = read_match_sum*h1_prob
	for k,sect_sum in s.items():
		if sect_sum > read_match_sum*bkg[k]:
			v2 = read_match_sum*bkg[k]*(1-bkg[k])
			m2 = read_match_sum*bkg[k]
			llr = np.log(v2**.5/v1**.5) + .5*((sect_sum-m2)**2/v2 - (sect_sum-m1)**2/v1)
			if llr > thresh:
				LLR.append((llr,k))
	return max(LLR)[1]

def get_blocks(Fc,fv,block_size):
	B = [bitarray() for _ in Fc]
	for b in range(len(B)):
		B[b].fromfile(Fc[b],block_size)
	V = np.fromfile(fv,dtype=np.float32,count=block_size*8)
	return B,V

def get_x(B,V,x):
	return [b for b in range(len(B)) if B[b][x]],V[x]

def best_partitions(cols,ids,Fcp,fvp,background_probs,hash_size):
	BestPartitions = []
	col_sort = np.argsort(cols)
	if len(col_sort) > 0:
		Fc = [open(fcp) for fcp in Fcp]
		fv = open(fvp)
		read_ids = np.empty(10**9,dtype=np.uint32)
		# only 16k clusters allowed here!
		cluster_ids = np.empty(10**9,dtype=np.int16)
		cluster_values = np.empty(10**9,dtype=np.float32)
		ix = 0
		last_x = None
		last_block = 0
		for i in col_sort:
			r_col = cols[i]
			r_id = ids[i]
			if r_col >= last_block:
				B,V = get_blocks(Fc,fv,block_size=min(5*10**9/len(Fc),(2**hash_size-last_block)/8))
				offset = last_block
				last_block += len(V)
			if r_col != last_x:
				col_clusters,col_value = get_x(B,V,int(r_col-offset))
				last_x = r_col
			if len(col_clusters) > 0:
				read_ids[ix] = r_id
				cluster_ids[ix] = -1
				cluster_values[ix] = col_value
				ix += 1
				for cluster in col_clusters:
					read_ids[ix] = r_id
					cluster_ids[ix] = cluster
					cluster_values[ix] = col_value
					ix += 1
		col_sort = None
		read_ids = read_ids[:ix]
		cluster_ids = cluster_ids[:ix]
		cluster_values = cluster_values[:ix]
		Fc = [fc.close() for fc in Fc]
		fv.close()
		id_sort = np.argsort(read_ids)
		if len(id_sort) > 0:
			current_id = read_ids[id_sort[0]]
			scores = defaultdict(float)
			for i in id_sort:
				if read_ids[i] != current_id:
					BestPartitions.append((current_id,max_log_lik_ratio(scores,background_probs)))
					current_id = read_ids[i]
					scores = defaultdict(float)
				scores[cluster_ids[i]] += cluster_values[i]
			BestPartitions.append((current_id,max_log_lik_ratio(scores,background_probs)))
	return BestPartitions

help_message = 'usage example: python write_partition_parts.py -r 1 -i /project/home/hashed_reads/ -o /project/home/cluster_vectors/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hr:i:o:',["--filerank=","inputdir=","outputdir="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-r',"--filerank"):
			fr = int(arg)-1
		elif opt in ('-i','--inputdir'):
			inputdir = arg
			if inputdir[-1] != '/':
				inputdir += '/'
		elif opt in ('-o','--outputdir'):
			outputdir = arg
			if outputdir[-1] != '/':
				outputdir += '/'
	hashobject = Fastq_Reader(inputdir,outputdir)
	cp = np.load(hashobject.output_path+'cluster_probs.npy')
	cluster_probs = dict(enumerate(cp))
	Cluster_Files = glob.glob(os.path.join(hashobject.output_path,'*.bits'))
	Cluster_Files = [(int(fp[fp.rfind('/')+1:fp.index('.bits')]),fp) for fp in Cluster_Files]
	Cluster_Files.sort()
	Cluster_Files = [fp[1] for fp in Cluster_Files]
	Hashq_Files = glob.glob(os.path.join(hashobject.input_path,'*.hashq.*'))
	Hashq_Files = [fp for fp in Hashq_Files if '.tmp' not in fp]
	Hashq_Files.sort()
	infile = Hashq_Files[fr]
	outpart = infile[-5:-3]
	sample_id = infile[infile.rfind('/')+1:infile.index('.hashq')]
	outfile = hashobject.output_path + sample_id + '.fastq.' + outpart + '.tmp'
	g_out_tmp = open(outfile,'w')
	f = gzip.open(infile)
	r_cols = np.empty(10**9,dtype=np.uint64)
	r_ids = np.empty(10**9,dtype=np.uint32)
	r_id = 0
	rx = 0
	for a in hashobject.hash_read_generator(f):
		if r_id%10**6 == 0:
			IdClusters = best_partitions(r_cols[:rx],r_ids[:rx],Cluster_Files,hashobject.output_path+'cluster_vals.npy',cluster_probs,hashobject.hash_size)
			for k,v in IdClusters:
				if v != None:
					g_out_tmp.write('%d\t%d\n' % (k,v))
			r_cols = np.empty(10**9,dtype=np.uint64)
			r_ids = np.empty(10**9,dtype=np.uint32)
			rx = 0
		for x in a[2]:
			r_cols[rx] = x
			r_ids[rx] = r_id
			rx += 1
		r_id += 1
	IdClusters = best_partitions(r_cols[:rx],r_ids[:rx],Cluster_Files,hashobject.output_path+'cluster_vals.npy',cluster_probs,hashobject.hash_size)
	for k,v in IdClusters:
		if v != None:
			g_out_tmp.write('%d\t%d\n' % (k,v))
	g_out_tmp.close()
	# pass over hashq again and write full reads to partition
	g_out_tmp = open(outfile)
	mapped_id,mapped_cluster = np.fromstring(g_out_tmp.readline(),dtype=np.uint64,sep='\t')
	f.seek(0)
	F = {}
	r_id = 0
	for a in hashobject.hash_read_generator(f):
		while mapped_id < r_id:
			try:
				mapped_id,mapped_cluster = np.fromstring(g_out_tmp.readline(),dtype=np.int64,sep='\t')
			except:
				mapped_id,mapped_cluster = (r_id+.5,None)
		if mapped_id == r_id:
			if mapped_cluster not in F:
				F[mapped_cluster] = open(hashobject.output_path+'.'.join([str(mapped_cluster),sample_id,'fastq',outpart]),'w')
			F[mapped_cluster].write(a[0]+'\n')
		r_id += 1
	f.close()
	g_out_tmp.close()
	os.system('rm '+outfile)
	for f in F.values():
		f.close()
	
