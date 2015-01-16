from random import randint
import glob,os
import numpy as np
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

	# A VERY HACKED METHOD FOR DETERMINING READ PAIR ID NAMEOLOGY (ie readid/1,readid/2 vs readid 1,readid 2)
	def id_type(self,f):
		initial_pos = f.tell()
		L = [f.readline() for _ in range(500)]
		L = [l for l in L if l]
		Ids = [l.strip().split() for l in L if l[0]=='@']
		pair_type = None
		for i in range(len(Ids)-1):
			if (Ids[i][0] == Ids[i+1][0]) and (len(Ids[i]) > 1) and (len(Ids[i+1]) > 1):
				if (Ids[i][1][0] == '1') and (Ids[i+1][1][0] == '2') and (Ids[i][1][1:] == Ids[i+1][1][1:]):
					pair_type = 1
					break
			elif (Ids[i][0][:-1] == Ids[i+1][0][:-1]):
				if ((Ids[i][0][-2:] == '/1') and (Ids[i+1][0][-2:] == '/2')) or ((Ids[i][0][-2:] == '.1') and (Ids[i+1][0][-2:] == '.2')):
					pair_type = 2
					break
		f.seek(initial_pos)
		return pair_type

	# Replace with hashq_read_mapper style that doesn't use read_until_new
	def hash_read_generator(self,file_object,max_reads=10**15,newline='\n'):
		line = file_object.readline().strip()
		lastlinechar = ''
		read_strings = []
		r = 0
		while (line != '') and (r < max_reads):
			# read_id (always) and quality (sometimes) begin with '@', but quality preceded by '+' 
			if (line[0] == '@') and (lastlinechar != '+'):
				if len(read_strings) == 5:
					try:
						I = newline.join(read_strings[:-1])
						B = np.fromstring(read_strings[-1][10:-2],dtype=np.uint64,sep=',')
						yield (I,B[0],B[1:])
					except Exception,err:
						print str(err)
					r += 1
				read_strings = []
			read_strings.append(line)
			lastlinechar = line[0]
			line = file_object.readline().strip()

	def read_generator(self,file_object,max_reads=10**15,verbose_ids=False,raw_reads=False):
		self.set_quality_codes(file_object)
		line = 'dummyline'
		r = 0
		while (line != '') and (r < max_reads):
			line = file_object.readline()
			if line:
				if line[0] == '@':
					try:
						# ASSUMING READ PAIRS ARE SPLIT INTO THEIR OWN LINES
						verbose_id = line
						I = line.strip()
						line = file_object.readline()
						verbose_id += line
						S = line.strip()
						# make sure this looks like sequence data - will break if upper and lower AaCcTtGg are used
						if len(set(S)) > 5:
							if len(set(S.upper())) > 5:
								raise Exception
							else:
								S = S.upper()
						verbose_id += file_object.readline()
						line = file_object.readline()
						verbose_id += line
						if raw_reads:
							yield verbose_id
						else:
							if verbose_ids:
								I = verbose_id
							Q = [self.quality_codes[c] for c in line.strip()]
							if (S) and (Q):
								low_qual = 0
								i = -1
								while (i < len(S)-self.kmer_size) and (low_qual < 3):
									if Q[i+self.kmer_size] < 3:
										low_qual += 1
									i += 1
								yield {'_id': I,'s': S[:i+self.kmer_size],'q': Q[:i+self.kmer_size]}
						r += 1
					except Exception,err:
						#print 'warning: fastq read_generator error',str(err)
						pass

	def set_quality_codes(self,f):
		last = f.tell()
		L = [f.readline() for _ in range(4000)]
		f.seek(last)
		L = [l for l in L if l]
		x = [sum([1 for l in L[i::4] if l[0]=='+']) for i in range(4)]
		x = x.index(max(x))+1
		o33 = 0
		o64 = 0
		for l in L[x::4]:
			for c in l.strip():
				oc = ord(c)
				if oc < 74:
					o33 += 1
				else:
					o64 += 1
		if 3*o33 > o64:
			self.quality_codes = dict([(chr(x),x-33) for x in range(33,33+94)])
		else:
			self.quality_codes = dict([(chr(x),x-64) for x in range(64-5,64+63)])

	def rand_kmers_for_wheel(self,total_kmers):
		RP = glob.glob(os.path.join(self.input_path,'*.fastq.*'))
		if len(RP) > 100:
			import random
			RP = random.sample(RP,100)
		elif len(RP) == 0:
			# single file per sample
			RP = glob.glob(os.path.join(self.input_path,'*.fastq'))
		kmers_per_file = max(total_kmers/len(RP),5)
		g = open(self.input_path+'random_kmers.fastq','w')
		total = 0
		for rp in RP:
			f = open(rp)
			kf = 0
			fails = 0
			while kf < kmers_per_file:
				try:
					g.write(self.rand_kmer(f))
					kf += 1
				except Exception,err:
					fails += 1
					print str(err)
					if fails > 100:
						break
			f.close()
			total += kf
			if total > total_kmers:
				break
		g.close()

	def rand_kmer(self,f,max_seek=10**8):
		while True:
			f.seek(randint(0,max_seek))
			# ultra slow - partly due to setting qual scores every time
			rs = [_ for _ in self.read_generator(f,max_reads=1,verbose_ids=True)]
			if len(rs) > 0:
				if len(rs[0]['s']) > self.kmer_size:
					break
			else:
				max_seek /= 10
			if max_seek == 0:
				raise Exception
		ri = min(20,randint(0,len(rs[0]['s'])-self.kmer_size))
		rs = rs[0]['_id'].split('\n')
		return '\n'.join([rs[0],rs[1][ri:ri+self.kmer_size],rs[2],rs[3][ri:ri+self.kmer_size]+'\n']) 
