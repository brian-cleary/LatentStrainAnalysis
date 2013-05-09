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
		L = [f.readline() for _ in range(15)]
		Ids = [l.strip().split() for l in L if l[0]=='@']
		pair_type = None
		for i in range(len(Ids)-1):
			if (Ids[i][0] == Ids[i+1][0]):
				if (Ids[i][1][0] == '1') and (Ids[i+1][1][0] == '2') and (Ids[i][1][1:] == Ids[i+1][1][1:]):
					pair_type = 1
					break
			elif (Ids[i][0][:-1] == Ids[i+1][0][:-1]):
				if (Ids[i][0][-2:] == '/1') and (Ids[i+1][0][-2:] == '/2'):
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

	def read_generator(self,file_object,max_reads=10**15,verbose_ids=False):
		self.set_quality_codes()
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
						line = file_object.readline()
						verbose_id += line
						S = line.strip()
						# make sure this looks like sequence data - will break if upper and lower AaCcTtGg are used
						if len(set(S)) > 5:
							raise Exception
						verbose_id += file_object.readline()
						line = file_object.readline()
						verbose_id += line
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

	def set_quality_codes(self):
		self.quality_codes = {'$': 3, '(': 7, ',': 11, '0': 15, '4': 19, '8': 23, '<': 27, '@': 31, 'D': 35, 'H': 39, 'L': 43, 'P': 47, 'T': 51, 'X': 55, '\\': 59, '`': 63, 'd': 67, 'h': 71, 'l': 75, 'p': 79, 't': 83, 'x': 87, '|': 91, '#': 2, "'": 6, '+': 10, '/': 14, '3': 18, '7': 22, ';': 26, '?': 30, 'C': 34, 'G': 38, 'K': 42, 'O': 46, 'S': 50, 'W': 54, '[': 58, '_': 62, 'c': 66, 'g': 70, 'k': 74, 'o': 78, 's': 82, 'w': 86, '{': 90, '"': 1, '&': 5, '*': 9, '.': 13, '2': 17, '6': 21, ':': 25, '>': 29, 'B': 33, 'F': 37, 'J': 41, 'N': 45, 'R': 49, 'V': 53, 'Z': 57, '^': 61, 'b': 65, 'f': 69, 'j': 73, 'n': 77, 'r': 81, 'v': 85, 'z': 89, '~': 93, '!': 0, '%': 4, ')': 8, '-': 12, '1': 16, '5': 20, '9': 24, '=': 28, 'A': 32, 'E': 36, 'I': 40, 'M': 44, 'Q': 48, 'U': 52, 'Y': 56, ']': 60, 'a': 64, 'e': 68, 'i': 72, 'm': 76, 'q': 80, 'u': 84, 'y': 88, '}': 92}

	def rand_kmers_for_wheel(self,total_kmers):
		RP = glob.glob(os.path.join(self.input_path,'*.fastq.*'))
		kmers_per_file = max(total_kmers/len(RP),5)
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
			rs = [_ for _ in self.read_generator(f,max_reads=1,verbose_ids=True)]
			if len(rs) > 0:
				if len(rs[0]['s']) > self.kmer_size:
					break
		ri = min(20,randint(0,len(rs[0]['s'])-self.kmer_size))
		rs = rs[0]['_id'].split('\n')
		return '\n'.join([rs[0],rs[1][ri:ri+self.kmer_size],rs[2],rs[3][ri:ri+self.kmer_size]+'\n']) 
