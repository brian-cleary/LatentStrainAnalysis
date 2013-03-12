from random import randint
import glob,os
from LSA import LSA
from hyper_sequences import Hyper_Sequences
from hash_counting import Hash_Counting
from cluster_analysis import Cluster_Analysis

class Fastq_Reader(Cluster_Analysis,Hash_Counting,Hyper_Sequences,LSA):

	def __init__(self,indir,outdir,hash_size=999,new_hash=(None,None)):
		super(Fastq_Reader,self).__init__(indir,outdir)
		if new_hash == (None,None):
			try:
				self.get_wheels(spoke_limit=hash_size,wheel_limit=1)
				self.hash_size = self.Wheels[-1]['s'] + 1
				self.kmer_size = len(self.Wheels[0]['p'])
			except:
				self.Wheels = None
				self.hash_size = None
				self.kmer_size = None
			self.newline_proxy = 'NEWLINE'
		else:
			self.hash_size = new_hash[0]
			self.kmer_size = new_hash[1]

	# Replace with hashq_read_mapper style that doesn't use read_until_new
	def hash_read_generator(self,file_object,max_reads=10**15,newline='\n'):
		line = file_object.readline().strip()
		read_strings = []
		r = 0
		while (line != '') and (r < max_reads):
			if line[0] == '@':
				if len(read_strings) == 5:
					try:
						I = newline.join(read_strings[:-1])
						B = [int(c) for c in read_strings[-1][10:-2].split(',')]
						yield (I,B[0],B[1:])
					except Exception,err:
						print str(err)
					r += 1
				read_strings = []
			read_strings.append(line)
			line = file_object.readline().strip()

	def read_generator(self,file_object,max_reads=10**15,verbose_ids=False,return_kmers=True):
		line = 'dummyline'
		r = 0
		while (line != '') and (r < max_reads):
			line = file_object.readline()
			if line:
				if line[0] == '@':
					try:
						# ASSUMING READ PAIRS ARE SPLIT INTO THEIR OWN LINES
						verbose_id = line
						#I = line.split()[1]
						#I = I[I.index(':')+1:]
						#I = line[line.index(':')+1:].strip()
						I = line.strip()
						verbose_id += file_object.readline()
						S = verbose_id.split('\n')[-2]
						verbose_id += file_object.readline()
						verbose_id += file_object.readline()
						if verbose_ids:
							I = verbose_id
						Q = self.quality_code_to_int(verbose_id.split('\n')[-2])
						if return_kmers:
							for kmer in self.kmers_from_read(S,Q,self.kmer_size):
								if kmer['q'].count(2) > 2:
									break
								if verbose_ids:
									kmer['_id'] = I
								yield kmer
						elif (S) and (Q):
							yield {'_id': I,'s': S,'q': Q}
						r += 1
					except Exception,err:
						print 'warning: fastq read_generator error',str(err)

	def quality_code_to_int(self,code_string):
		quality_codes = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
		return [quality_codes.index(c) for c in code_string]

	def kmers_from_read(self,s,q,k):
		i = 0
		while i < len(s)-k+1:
			yield {'_id': i,'s': s[i:i+k],'q': q[i:i+k]}
			i += 1

	def rand_kmers_for_wheel(self,total_kmers):
		RP = glob.glob(os.path.join(self.input_path,'*.fastq.*'))
		kmers_per_file = total_kmers/len(RP)
		g = open(self.input_path+'random_kmers.fastq','w')
		for rp in RP:
			f = open(rp)
			kf = 0
			while kf < kmers_per_file:
				try:
					g.write(self.rand_kmer(f))
					kf += 1
				except Exception,err:
					print str(err)
			f.close()
		g.close()

	def rand_kmer(self,f,max_seek=10**8):
		while True:
			f.seek(randint(0,max_seek))
			rs = [_ for _ in self.read_generator(f,max_reads=1,verbose_ids=True,return_kmers=False)]
			if len(rs) > 0:
				if len(rs[0]['s']) > self.kmer_size:
					break
		ri = randint(0,len(rs[0]['s'])-self.kmer_size)
		rs = rs[0]['_id'].split('\n')
		return '\n'.join([rs[0],rs[1][ri:ri+self.kmer_size],rs[2],rs[3][ri:ri+self.kmer_size]+'\n']) 
