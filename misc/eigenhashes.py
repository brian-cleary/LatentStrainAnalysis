import numpy as np
from scipy.sparse.linalg import svds
from scipy.sparse import csr_matrix
from itertools import combinations,islice
from multiprocessing import Pool
from scipy.spatial import distance
from LSA import LSA
from hyper_sequences import Hyper_Sequences
from hash_counting import Hash_Counting

class Eigenhashes(Hash_Counting,Hyper_Sequences,LSA):

	def __init__(self,inpath,outpath):
		super(Eigenhashes,self).__init__(inpath,outpath)
		self.get_wheels(wheel_limit=1)
		self.hash_size = self.Wheels[-1]['s'] + 1
		self.kmer_size = len(self.Wheels[0]['p'])
		self.pool = Pool()

	def matrix_from_file_paths(self,path_list):
		M = np.zeros((len(path_list),2**self.hash_size),dtype=np.float32)
		for p in range(len(path_list)):
			M[p,:] = self.open_count_hash(path_list[p])
		return M

	def nonzero_abundances(self,M):
		NZ = np.array(self.nonzero_cols(M),dtype=np.uint32)
		np.save(self.output_path+'nonzero_indices.npy',NZ)
		return M[:,NZ],NZ

	def nonzero_cols(self,M):
		i = 0
		Nonzeros = []
		while i < M.shape[1]:
			m = M[:,i:i+1000000].sum(0)
			for s in range(len(m)):
				if m[s] != 0:
					Nonzeros.append(i+s)
			i += 1000000
		return Nonzeros

	def calculate_global_weights(self,M,NZ):
		GWnz = self.global_entropy_weights(M)
		np.save(self.output_path+'global_nonzero_weights.npy',GWnz)
		GW = np.zeros(2**self.hash_size,dtype=np.float32)
		for i in xrange(len(NZ)):
			GW[NZ[i]] = GWnz[i]
		NZ = None
		np.save(self.output_path+'global_weights.npy',GW)
		GW = None
		return GWnz

	def conditioned_nonzeros(self,M,GWnz):
		M = np.log(M + 1)*GWnz
		GWnz = None
		np.save(self.output_path+'conditioned_nonzeros.npy',M)
		return M

	def global_entropy_weights(self,M,block_size=1000000):
		rows,cols = M.shape
		log_rows = np.log(rows)
		Global_Weights = np.zeros(cols,dtype=np.float32)
		block = []
		for ci,c in enumerate(M.T):
			block.append((ci,c,log_rows))
			if len(block) > block_size:
				updates = self.pool.map(pooled_entropy,block,max(1,len(block)/20))
				for u in updates:
					if u:
						Global_Weights[u[0]] = u[1]
				block = []
				updates = None
		if block:
			updates = self.pool.map(pooled_entropy,block,max(1,len(block)/20))
			for u in updates:
				if u:
					Global_Weights[u[0]] = u[1]
			block = []
			updates = None
		self.pool.close()
		self.pool.join()
		return Global_Weights

	# Note: explicitly creating M as a sparse matrix can consume massive memory as M becomes dense. However, svds (as opposed to np.linalg.svd) can still be used as a memory efficient, but slow, decomposition method.
	def eigenkmers(self,M,num_dims=10):
		L,V,R = svds(M,num_dims)
		np.save(self.output_path+'eigenleft.npy',L)
		np.save(self.output_path+'eigenvalues.npy',V)
		np.save(self.output_path+'eigenright.npy',R)
		return np.transpose(np.transpose(R)*V)

	def random_cols(self,M,n):
		num_cols = M.shape[1]
		RC = []
		while len(RC) < n:
			r = np.random.randint(0,num_cols)
			if abs(sum(M[:,r])) > 10**-9:
				RC.append(M[:,r])
		return RC

	def kmer_clusters(self,M,initial_clusters=200,cluster_thresh=0.8,cluster_iters=3,block_size=100):
		Clusters = []
		Centers = []
		for v in self.random_cols(M,initial_clusters):
			Clusters.append([])
			if len(Centers) > 0:
				Centers = np.concatenate((Centers,[v]))
			else:
				Centers = [v]
		num_cols = M.shape[1]
		for _ in range(cluster_iters):
			Clusters,Centers = self.cluster_centers(Clusters,Centers,M)
			if _ == cluster_iters-1:
				maxfits = 5
			else:
				maxfits = 1
			i = 0
			while i < num_cols:
				Clusters,Centers = self.distance_block((i,i+block_size),M,Clusters,Centers,cluster_thresh,max_cluster_fits=maxfits)
				i += block_size
				if i%10000000 == 0:
					print _,i,len(Clusters)
		return Clusters

	def distance_block(self,indices,M,C,Cm,ct,max_cluster_fits=1):
		Msum = M[:,indices[0]:indices[1]].sum(0)
		nonzero_indices = [indices[0]+i for i in range(len(Msum)) if abs(Msum[i])>10**-9]
		D = distance.cdist(np.transpose(M[:,nonzero_indices]),Cm,'cosine')
		for i in xrange(D.shape[0]):
			found_cluster = False
			MI = D[i,:].argsort()[:max_cluster_fits]
			for min_i in MI:
				if D[i,min_i] < 1-ct:
					C[min_i].append(nonzero_indices[i])
					found_cluster = True
			if not found_cluster:
				C.append([nonzero_indices[i]])
				Cm = np.concatenate((Cm,[M[:,nonzero_indices[i]]]))
		return C,Cm

	def cluster_centers(self,C,Cm,M,combine_thresh=.8):
		for k in range(len(C)):
			v = C[k]
			if len(v) > 0:
				sampled_members = np.array(self.random_cols(M[:,v],min(150000,len(v))))
				Cm[k,:] = sampled_members.sum(axis=0)/len(sampled_members)
		remove_clusters = {}
		D = distance.pdist(Cm,'cosine')
		D = D < (1 - combine_thresh)
		i = 0
		j = 1
		for d in D:
			if j >= Cm.shape[0]:
				i += 1
				j = i+1
			if d:
				if len(C[i]) >= len(C[j]):
					remove_clusters[j] = True
				else:
					remove_clusters[i] = True
			j += 1
		return [[] for _ in range(len(C)-len(remove_clusters))],Cm[[i for i in range(len(C)) if i not in remove_clusters],:]

	def save_clusters(self,Clusters,Nonzeros):
		for i in range(len(Clusters)):
			c = [Nonzeros[x] for x in Clusters[i]]
			np.save(self.output_path+str(i)+'.cluster.npy',c)
		np.save(self.output_path+'kmer_cluster_sizes.npy',[len(c) for c in Clusters])

def pooled_entropy(args):
	col_index,col,log_rows = args
	gf = sum(col)
	rc_updates = []
	if gf > 0:
		g_weight = 1 + sum([tf/gf*np.log(tf/gf)/log_rows for tf in col if tf>0])
		return (col_index,g_weight)
