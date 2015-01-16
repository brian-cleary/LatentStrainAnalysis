import numpy as np
from scipy.spatial import distance
from gensim import models
import logging
from collections import defaultdict
from multiprocessing import Pool
from operator import itemgetter
from LSA import LSA
from hyper_sequences import Hyper_Sequences
from hash_counting import Hash_Counting

class StreamingEigenhashes(Hash_Counting,Hyper_Sequences,LSA):

	def __init__(self,inpath,outpath,get_pool=8):
		super(StreamingEigenhashes,self).__init__(inpath,outpath)
		self.get_wheels(wheel_limit=1)
		self.hash_size = self.Wheels[-1]['s'] + 1
		self.kmer_size = len(self.Wheels[0]['p'])
		if get_pool > 0:
			self.pool = Pool(get_pool)

	def corpus_idf_from_hash_paths(self,path_list):
		class NewCorpus(object):
			def __iter__(newself):
				for fp in path_list:
					nonzeros = np.load(fp)
					yield nonzeros
		return NewCorpus()

	def train_tfidf(self,corpus):
		doc_freqs = np.ones(2**self.hash_size,dtype=np.float32)
		total_docs = 2.
		for doc in corpus:
			doc_freqs[doc] += 1
			total_docs += 1.
		self.global_weights = np.log2(total_docs/doc_freqs)

	def kmer_corpus_to_disk(self,fp,mask=[]):
		H = np.array(self.open_count_hash(fp),dtype=np.float32)
		H[mask] = 0
		del mask
		norm = np.linalg.norm(H)/len(H)**.5
		X = np.memmap(fp+'.conditioned',dtype=np.float32,mode='w+',shape=(len(H),))
		self.global_weights = np.load(self.output_path+'global_weights.npy')
		X[:] = H*self.global_weights/norm
		del X

	def kmer_corpus_from_disk(self,o=None):
		if o == None:
			o = (0,2**self.hash_size)
		o = (o[0],min(o[1],2**self.hash_size - o[0]))
		F = [open(self.path_dict[i],'rb') for i in range(len(self.path_dict))]
		for f in F:
			f.seek(o[0]*4)
		class NewCorpus(object):
			def __iter__(newself):
				rows_processed = 0
				while True:
					H = self.new_chunk(F,n=min(10000,o[1]-rows_processed))
					for h in H:
						nz = np.nonzero(h)[0]
						yield [(x,h[x]) for x in nz]
					rows_processed += H.shape[0]
					if (rows_processed >= o[1]) or (H.shape[0] == 0):
						break
				for f in F:
					f.close()
		return NewCorpus()

	def new_chunk(self,F,n=10000):
		n = int(n)
		X = np.zeros((n,len(F)),dtype=np.float32)
		for j in range(len(F)):
			X[:,j] = np.fromfile(F[j],dtype=np.float32,count=n)
		return X

	def train_kmer_lsi(self,kmer_corpus,num_dims=20,single=False):
		logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
		if single:
			return models.LsiModel(kmer_corpus,num_topics=num_dims,id2word=self.path_dict,distributed=False,chunksize=200000)
		else:
			return models.LsiModel(kmer_corpus,num_topics=num_dims,id2word=self.path_dict,distributed=True,chunksize=200000)

	def lsi_cluster_index(self,lsi,random_chunk=0.00002,cluster_iters=200):
		Clusters = {}
		Index = np.zeros((0,lsi.num_topics))
		chunk_size = random_chunk*2**self.hash_size
		for ci in range(cluster_iters):
			seed_vectors = lsi[self.kmer_corpus_from_disk(o=(np.random.randint(0,2**self.hash_size-chunk_size),chunk_size))]
			Clusters,Index = self.merge_index(seed_vectors,Index,Clusters)
			Clusters,Index = self.collapse_index(Index,Clusters)
			print ci,len(Clusters)
		return Index

	def lsi_cluster_part(self,offsets,lsi,Index):
		cluster_thresh = self.cluster_thresh
		Clusters = [np.empty((offsets[1],),dtype=np.int64) for _ in range(Index.shape[0])]
		Sizes = np.zeros(Index.shape[0],dtype=np.int64)
		num_best = 5
		all_vectors = lsi[self.kmer_corpus_from_disk(o=offsets)]
		vector_block = []
		block_index = []
		for col,doc in enumerate(all_vectors):
			# this is slow and not particularly clever...
			if len(doc) > 0:
				block_index.append(col+offsets[0])
				a = np.zeros(Index.shape[1])
				for x in doc:
					a[x[0]] = x[1]
				vector_block.append(a)
			if len(vector_block) == 10**4:
				D = distance.cdist(vector_block,Index,'cosine')
				for indexed_doc in enumerate(D):
					colx,fits = indexed_doc
					MI = fits.argsort()[:num_best]
					for clust in MI:
						if fits[clust] < 1-cluster_thresh:
							Clusters[clust][Sizes[clust]] = block_index[colx]
							Sizes[clust] += 1
						else:
							break
				block_index = []
				vector_block = []
		if len(vector_block) > 0:
			D = distance.cdist(vector_block,Index,'cosine')
			for indexed_doc in enumerate(D):
				colx,fits = indexed_doc
				MI = fits.argsort()[:num_best]
				for clust in MI:
					if fits[clust] < 1-cluster_thresh:
						Clusters[clust][Sizes[clust]] = block_index[colx]
						Sizes[clust] += 1
					else:
						break
		return [c[:Sizes[i]] for i,c in enumerate(Clusters)]

	def lsi_kmer_clusters(self,lsi,Index,c0,cs=10**5):
		cluster_thresh = self.cluster_thresh
		Clusters = [np.empty((cs,),dtype=np.int64) for _ in range(Index.shape[0])]
		Sizes = np.zeros(Index.shape[0],dtype=np.int64)
		#for j in range(0,2**self.hash_size,10**7):
		block = [((i,10**5),lsi,Index,cluster_thresh,self.input_path,self.output_path,self.path_dict) for i in range(j,min(j+10**7,2**self.hash_size),10**5)]
		results = self.pool.map(distance_pool,block)
		for r in results:
			for i in range(len(Clusters)):
				for x in enumerate(r[i]):
					Clusters[i][Sizes[i]+x[0]] = x[1]
				Sizes[i] += len(r[i])
		self.pool.close()
		self.pool.join()
		return [Clusters[i][:Sizes[i]] for i in range(len(Clusters))]

	def merge_index(self,V,I,C):
		thresh = self.cluster_thresh
		MergeFits = defaultdict(list)
		for doc in V:
			if len(doc) > 0:
				a = np.zeros(I.shape[1])
				for x in doc:
					a[x[0]] = x[1]
				if I.shape[0] > 0:
					fits = distance.cdist([a],I,'cosine')[0]
					clust = fits.argsort()[0]
					if fits[clust] < 1-thresh:
						MergeFits[clust].append(a)
					else:
						I = np.concatenate((I,[a]))
						C[len(C)] = 1
				else:
					I = np.array([a])
					C[0] = 1
		for k,v in MergeFits.items():
			I[k,:] = np.concatenate(([I[k,:]*C[k]],v)).sum(0)/(len(v)+C[k])
			C[k] += len(v)
		return C,I

	def collapse_index(self,I,C):
		combine_thresh = self.cluster_thresh
		remove_clusters = {}
		D = distance.pdist(I,'cosine')
		D = D < (1 - combine_thresh)
		i = 0
		j = 1
		for d in D:
			if j >= I.shape[0]:
				i += 1
				j = i+1
			if d:
				if C[i] >= C[j]:
					remove_clusters[j] = True
				else:
					remove_clusters[i] = True
			j += 1
		Cnew = {}
		for i in range(len(C)):
			if i not in remove_clusters:
				Cnew[len(Cnew)] = C[i]
		return Cnew,I[[i for i in range(len(C)) if i not in remove_clusters],:]