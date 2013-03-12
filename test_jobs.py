#!/usr/bin/env python

import numpy as np
from scipy.sparse.linalg import svds
from scipy.sparse import csr_matrix
from itertools import combinations
from multiprocessing import Pool
from scipy.spatial import distance
from bitarray import bitarray
from ctypes import c_uint8
from time import time
import itertools

def doit(args):
	x,y = args
	z = sum([np.random.randn(1) for _ in range(1000)])
	return x*y*z

if __name__ == "__main__":
	pool = Pool()
	M = np.zeros((20,2**29),dtype=np.float32)
	M[np.random.randint(0,M.shape[0]),np.random.randint(0,M.shape[1])] = 8
	#M = csr_matrix(M)
	print M
	t0 = time()
	results = []
	block = []
	for ci,c in enumerate(M[:,:100000].T):
		block.append((ci,c))
		if len(block) > 10000:
			updates = pool.map(doit,block,10000/20)
			if updates:
				results.extend(updates)
			block = []
	pool.close()
	pool.join()
	print 'total time:',time()-t0
	print 'total results:',len(results)