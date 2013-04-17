from bitarray import bitarray
from ctypes import c_uint16
import glob,os
from collections import defaultdict
import numpy as np
import scipy.stats as stats
import gzip
from LSA import LSA

class Hash_Counting(LSA):

	def __init__(self,inputpath,outputpath):
		super(Hash_Counting,self).__init__(inputpath,outputpath)

	def hash_counts_from_hashq(self,fileprefix,multi_files=False):
		H = (c_uint16*2**self.hash_size)()
		if multi_files:
			FP = glob.glob(os.path.join(self.output_path,fileprefix+'.*.hashq.*'))
		else:
			FP = [self.output_path+fileprefix]
		for filename in FP:
			f = gzip.open(filename)
			last = None
			while last != f.tell():
				last = f.tell()
				for a in self.hash_read_generator(f):
					for b in a[2]:
						H[b] = min(65535,H[b]+1)
			f.close()
		f0 = open(self.output_path+fileprefix+'.count.hash','wb')
		f0.write(H)
		f0.close()
		return H

	def open_count_hash(self,file_path):
		f = open(file_path,'rb')
		H = (c_uint16*2**self.hash_size)()
		f.readinto(H)
		f.close()
		return H

	def bitarray_from_array(self,A):
		H = bitarray(2**self.hash_size)
		H.setall(False)
		for a in A:
			H[a] = True
		return H

	def membership_generator(self,H,Hkeys,outsuffix,match_thresh=3.84,h1_prob=0.8):
		f = gzip.open(self.infile)
		g = gzip.open(self.outfile + '.' + str(outsuffix),'wb')
		for a in self.hash_read_generator(f,newline=self.newline_proxy):
			try:
				read_set = set(a[2])
				read_match_sum = self.global_weights[a[2]].sum(dtype=np.float64)
				for h in range(len(H)):
					sect_sum = self.global_weights[list(read_set & H[h][0])].sum(dtype=np.float64)
					if sect_sum > read_match_sum*H[h][1]:
						# we could also round to integers, then use binom.pmf
						# log ratio of Pr( instersection size | cluster size, p=h1_prob ) to Pr( intersection size | cluster size, p=(cluster size)/(total size))
						log_likelihood_ratio = np.log(stats.norm.pdf(sect_sum,read_match_sum*h1_prob,(read_match_sum*h1_prob*(1-h1_prob))**.5)/stats.norm.pdf(sect_sum,read_match_sum*H[h][1],(read_match_sum*H[h][1]*(1-H[h][1]))**.5))
						if log_likelihood_ratio > match_thresh:
							g.write('%s\t%f\t%s' % (Hkeys[h],log_likelihood_ratio,a[0]))
							g.write('\n')
			except:
				pass
		g.close()
		f.close()
	
	# with many clusters this may throw an error for too many open files
	def fastq_from_intermediate_output(self,group,sample_id):
		PF = glob.glob(os.path.join(self.input_path,group+'.*'))
		Reads = {}
		F = {}
		for pf in PF:
			f = gzip.open(pf)
			for l in f:
				try:
					l_id,l_score,l_info = l.strip().split('\t')
					read_id = l_info.split(self.newline_proxy)[0]
					if read_id not in Reads:
						Reads[read_id] = (l_info,[])
					Reads[read_id][1].append((float(l_score),l_id))
				except Exception,err:
					pass
			f.close()
		for read_id,values in Reads.iteritems():
			read_info,partitions = values
			top_partition = max(partitions)[1]
			if top_partition not in F:
				F[top_partition] = open(self.output_path+'.'.join([top_partition,sample_id,'fastq',group]),'w')
			F[top_partition].write(read_info.replace(self.newline_proxy,'\n') + '\n')
		for f in F.values():
			f.close()

	def collision_report(self):
		f = gzip.open(self.infile)
		B = defaultdict(list)
		N = 0
		for a in self.hash_read_generator(f,max_reads=10**5):
			try:
				seq = a[0].split('\n')[1]
				for b in range(len(a[2])):
					B[a[2][b]].append(seq[b:b+self.kmer_size])
					N += 1
			except:
				pass
		f.close()
		D = []
		for v in B.values():
			for i in range(len(v)):
				for j in range(i+1,len(v)):
					D.append(sum(a!=b for a,b in zip(v[i],v[j])))
		return (N,len(D),np.histogram(D,bins=np.arange(self.kmer_size)))
