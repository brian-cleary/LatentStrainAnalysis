import numpy as np
import cPickle
from operator import itemgetter
from LSA import LSA

class Hyper_Sequences(LSA):

	def __init__(self,inputpath,outputpath):
		super(Hyper_Sequences,self).__init__(inputpath,outputpath)

	def generator_to_coords(self,sequence_generator):
		for s in sequence_generator:
			coords = self.letters_to_coords(s)
			yield s['_id'],coords

	def letters_to_coords(self,S):
		# THIS EQUATES A/T AND C/G MISMATCHES WITH NEGATIVE MATCHES, AND ALL OTHERS AS NON-MATCHES
		ltc_q = {'A': complex(-1,0),'T': complex(1,0),'C': complex(0,-1),'G': complex(0,1)}
		p_correct = 1-1./2000
		ltc_no_q = {'A': complex(-1,0)*p_correct,'T': complex(1,0)*p_correct,'C': complex(0,-1)*p_correct,'G': complex(0,1)*p_correct}
		if 'q' in S:
			return np.array([ltc_q.get(l,complex(0,0)) for l in S['s']])*self.quality_to_prob(S['q'])
		else:
			return [ltc_no_q.get(l,complex(0,0)) for l in S['s']]

	def quality_to_prob(self,Q):
		return np.array([1-10**(-q/10.) for q in Q])

	def set_wheels(self,wheels=200):
		random_kmer_path = self.input_path + 'random_kmers.fastq'
		Wheels = []
		for w in xrange(wheels):
			Wheels += self.one_wheel(w,random_kmer_path)
		Wheels.sort()
		f = open(self.output_path+'Wheels.txt','w')
		cPickle.dump(Wheels,f)
		f.close()

	def get_wheels(self,spoke_limit=999,wheel_limit=999999):
		try:
			f = open(self.output_path+'Wheels.txt')
		except:
			f = open(self.input_path+'Wheels.txt')
		Wheels = cPickle.load(f)
		f.close()
		self.Wheels = [{'w': x[0],'s': x[1],'p': x[2],'c': x[3]} for x in Wheels if (x[0] < wheel_limit) and (x[1] < spoke_limit)]

	def coords_to_bins(self,A,C,reverse_compliments=True):
		num_wheels = self.Wheels[-1]['w'] + 1
		num_spokes = self.Wheels[-1]['s'] + 1
		pow2 = np.array([2**j for j in range(num_spokes)])
		Wc = np.array([w['c'] for w in self.Wheels])
		C = np.array(C)
		L = np.dot(C,np.transpose([w['p'] for w in self.Wheels]).conjugate())
		B = [np.dot((L[:,ws:ws+num_spokes] > Wc[ws:ws+num_spokes]),pow2) for ws in range(0,num_wheels*num_spokes,num_spokes)]
		if reverse_compliments:
			L = np.dot(C[:,::-1]*-1,np.transpose([w['p'] for w in self.Wheels]).conjugate())
			B2 = [np.dot((L[:,ws:ws+num_spokes] > Wc[ws:ws+num_spokes]),pow2) for ws in range(0,num_wheels*num_spokes,num_spokes)]
			return A,self.pick_one_from_rc_pair(B,B2)
		else:
			return A,B

	# ONLY NEED TO WRITE DOWN ONE KMER FROM REVERSE COMPLIMENT PAIR
	def pick_one_from_rc_pair(self,b1,b2,mx=1000000):
		B = [np.array([b1[i],b2[i]]) for i in range(len(b1))]
		return [b[np.mod(b,mx).argmin(0),range(b.shape[1])] for b in B]

	def generator_to_bins(self,sequence_generator,rc=True):
		C = []
		A = []
		for a,c in self.generator_to_coords(sequence_generator):
			for i in range(len(c)-self.kmer_size+1):
				A.append(a)
				C.append(c[i:i+self.kmer_size])
		if len(A) > 0:
			return self.coords_to_bins(A,C,reverse_compliments=rc)
		else:
			return None,None

	def one_wheel(self,w,rp):
		S = []
		f = open(rp)
		for s in range(self.hash_size):
			L = self.pick_leaf_noloc(self.kmer_size,f)
			P = self.affine_hull(L.values())
			C = P.pop()
			S.append((w,s,P,C))
		f.close()
		return S

	def pick_leaf_noloc(self,nodes,f):
		new_leaf = {}
		nl = [_ for _ in self.generator_to_coords(self.read_generator(f,max_reads=nodes))]
		for nlx in nl:
			new_leaf[len(new_leaf)] = list(nlx[1])
		return new_leaf

	def affine_hull(self,linear_system):
		# linear_system: d dimensions of n docs in this hyperplane
		for row in linear_system:
			row.append(-1)
		linear_system.append([0]*len(linear_system[0]))
		linear_system = np.array(linear_system)
		U,W,V = np.linalg.svd(linear_system)
		return list(V[-1,:])
