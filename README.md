Welcome to the Latent Strain Analysis (LSA) code repository!

LSA was developed as a pre-assembly tool for partitioning metagenomic reads. It uses a hyperplane hashing function and streaming SVD in order to find covariance relations between k-mers. The code, and the process outline in LSFScripts in particular, have been optimized to scale to massive data sets in fixed memory with a highly distributed computing environment.

Dependencies
	Python (2.7)
	NumPy
	SciPy
	Gensim
	GNU Parallel (http://www.gnu.org/software/parallel/) for the single-instance version
	Pyro4 for the LSF version

GETTING STARTED

To get started, download the repo (eg via git clone). Next, unpack testData.tar.gz (into original_reads/). This folder contains a subsample of 10k reads from each of the 18 metagenomic samples used in Sharon et. al. (SRA052203), and in the original LSA methods paper. Each sample has been randomly spiked with 0-2,000 mock reads from a Bacillus thuringiensis plasmid (NG_035027.1).

Begin by performing some hashing tasks using 6 threads, k-mers of length 33, and a hash size of 2^22:

	$ bash HashCounting.sh 6 33 22

Next, calculate the SVD and cluster the hashed k-mers using 6 threads, a hash size of 2^22 and a cluster threshold of 0.8:

	$ bash KmerSVDClustering.sh 6 22 .8

Finally, partition the original reads (4 threads):

	$ bash ReadPartitioning.sh 4

Since the process includes several points of randomization, the results from each run will vary. However, when the read partitioning step completes, a report is made to indicate the specificity and completeness of spiked read partitioning:

		Quantifying spike enrichment
		Total spiked reads: 21528
		Spiked read counts by partition (top 5)
		partition 5: 19771 spiked reads out of 20291 total reads in partition 5 (97%)
		partition 146: 1207 spiked reads out of 1557 total reads in partition 146 (77%)
		partition 131: 356 spiked reads out of 670 total reads in partition 131 (53%)
		partition 112: 142 spiked reads out of 492 total reads in partition 112 (28%)
		partition 129: 4 spiked reads out of 506 total reads in partition 129 (0%)

ANALYZING LARGER COLLECTIONS

LSA has been written to be highly efficient in analyzing very large collections. For data sets larger than a few Gb, significant gains in wall time can be made by running in a cluster environment. In these cases, the process is essentially the same as what is outlined above. Detailed steps and job array submission scripts can be found in LSFScripts/README.md.