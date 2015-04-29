Latent Strain Analysis Documentation
====================================

Welcome to the Latent Strain Analysis (LSA) documentation!

LSA was developed as a pre-assembly tool for partitioning metagenomic reads. It uses a hyperplane hashing function and streaming SVD in order to find covariance relations between k-mers. The code, and the process outlined in "Analyzing Large Data Sets" in particular, have been optimized to scale to massive data sets in fixed memory with a highly distributed computing environment.

Overview
--------

LSA operates at the level of k-mers, and, more specifically, on hashed k-mers. In the first steps of LSA we generate a hash function that maps k-mers to columns in a matrix. The rows of this matrix will represent different samples, and we are interested in the covariance structure of k-mers across samples, with the idea that this covariance can reflect the physical linkage between k-mers found in the same genome.

We can get at the covariance information via Singular Value Decomposition, and, since our matrices might be super huge to accomodate a large diversity of k-mers, we implement a streaming SVD so that the whole thing works in fixed memory. Then, working in the eigenspace of k-mer covariation, we find a k-mer clustering that is subsequently used to partition reads into disjoint subsets.

Contents
--------

.. toctree::
   :maxdepth: 2

   license
   getting_started
   large_collections


