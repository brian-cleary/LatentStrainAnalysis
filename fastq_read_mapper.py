#!/usr/bin/env python

import sys
import numpy as np
import cPickle

def read_generator(file_object,max_reads=10**15,verbose_ids=False,from_hashq=False,kmer_size=None):
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
					Q = quality_code_to_int(verbose_id.split('\n')[-2])
					if kmer_size:
						for kmer in kmers_from_read(S,Q,kmer_size):
							if kmer['q'].count(2) > 2:
								break
							#if verbose_ids:
							kmer['_id'] = I
							yield kmer
					elif from_hashq:
						B = RS[l]
						l += 1
						if 'k, bins:' in B:
							B = [int(c) for c in B[10:-2].split(',')]
							yield (I,B[0],B[1:])
					elif (S) and (Q):
						yield {'_id': I,'s': S,'q': Q}
					r += 1
				except:
					pass

quality_codes = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
def quality_code_to_int(code_string):
	return [quality_codes.index(c) for c in code_string]

def kmers_from_read(s,q,k):
	i = 0
	while i < len(s)-k+1:
		yield {'_id': i,'s': s[i:i+k],'q': q[i:i+k]}
		i += 1

def generator_to_bins(sequence_generator,Wheels,rc=False,return_terminals=False):
	C = []
	A = []
	for a,c in generator_to_coords(sequence_generator):
		if return_terminals:
			a = (a,c[0],c[-1])
		A.append(a)
		C.append(c)
	return coords_to_bins(A,C,Wheels,reverse_compliments=rc)

def generator_to_coords(sequence_generator):
	for s in sequence_generator:
		coords = letters_to_coords(s)
		yield s['_id'],coords

# THIS EQUATES A/T AND C/G MISMATCHES WITH NEGATIVE MATCHES, AND ALL OTHERS AS NON-MATCHES
ltc_q = {'A': complex(-1,0),'T': complex(1,0),'C': complex(0,-1),'G': complex(0,1)}
p_correct = 1-1./2000
ltc_no_q = {'A': complex(-1,0)*p_correct,'T': complex(1,0)*p_correct,'C': complex(0,-1)*p_correct,'G': complex(0,1)*p_correct}
def letters_to_coords(S):
	if 'q' in S:
		return np.array([ltc_q.get(l,complex(0,0)) for l in S['s']])*quality_to_prob(S['q'])
	else:
		return [ltc_no_q.get(l,complex(0,0)) for l in S['s']]

def quality_to_prob(Q):
	return np.array([1-10**(-q/10.) for q in Q])

def coords_to_bins(A,C,Wheels,reverse_compliments=False):
	if len(C) > 0:
		num_wheels = Wheels[-1]['w'] + 1
		num_spokes = Wheels[-1]['s'] + 1
		pow2 = [2**j for j in range(num_spokes)]
		L = np.dot(C,np.transpose([w['p'] for w in Wheels]).conjugate())
		L -= [w['c'] for w in Wheels]
		L = np.int_((np.sign(L) + 1)/2)
		B = [np.dot(L[:,ws:ws+num_spokes],pow2) for ws in range(0,num_wheels*num_spokes,num_spokes)]
		if reverse_compliments:
			L = np.dot(np.array(C)[:,::-1]*-1,np.transpose([w['p'] for w in Wheels]).conjugate())
			L -= [w['c'] for w in Wheels]
			L = np.int_((np.sign(L) + 1)/2)
			B2 = [np.dot(L[:,ws:ws+num_spokes],pow2) for ws in range(0,num_wheels*num_spokes,num_spokes)]
			return A,[[pick_one_from_rc_pair(B[i][j],B2[i][j]) for j in range(len(B[i]))] for i in range(len(B))]
		else:
			return A,B
	else:
		return None,[]

# ONLY NEED TO WRITE DOWN ONE KMER FROM REVERSE COMPLIMENT PAIR
def pick_one_from_rc_pair(b1,b2,mx=1000000):
	if (b1 % mx) < (b2 % mx):
		return b1
	else:
		return b2

def read_kmer_bins(k,W):
	A = []
	while A != None:
		try:
			A,B = generator_to_bins(read_generator(sys.stdin,max_reads=1000,kmer_size=k),W,rc=True)
			for b in range(len(B)):
				for a in range(len(A)):
					#print '%s\t%d' % (B[b][a],1)
					a_id = A[a].split()[0]
					print '%s\t%d' % (a_id,B[b][a])
					print '%s\t%d' % (mate_pair_id(a_id),B[b][a])
		except Exception,err:
			pass

def mate_pair_id(s):
	if s[-1] == '1':
		return s[:-1] + '2'
	else:
		return s[:-1] + '1'

def get_wheels(wheel_path,spoke_limit=999,wheel_limit=999999):
	f = open(wheel_path)
	Wheels = cPickle.load(f)
	f.close()
	Wheels = [{'w': x[0],'s': x[1],'p': x[2],'c': x[3]} for x in Wheels if (x[0] < wheel_limit) and (x[1] < spoke_limit)]
	return Wheels

if __name__ == "__main__":
	W = get_wheels('Wheels.txt',spoke_limit=28,wheel_limit=1)
	kmer_size = len(W[0]['p'])
	read_kmer_bins(kmer_size,W)
