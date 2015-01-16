Welcome to the Latent Strain Analysis (LSA) code repo!

LSA was developed in Python 2.7, and requires NumPy, SciPy, and Gensim (for streaming SVD calculation). Pyro4 is also required for running distributed, streaming SVD calculations, which are used in the LSF version of the tool. For the single instance version used in "Getting Started", GNU parallel (http://www.gnu.org/software/parallel/) is needed to make use of multiple cores.

GETTING STARTED

To get started, download the repo (eg via git clone). Nect, unpack testData.tar.gz (into original_reads/). This folder contains a subsample of 10k reads from each of the 18 metagenomic samples used in Sharon et. al. (SRA052203), and in the original LSA methods paper. Each sample has been randomly spiked with 0-2,000 mock reads from a Bacillus thuringiensis plasmid (NG_035027.1).

Begin by performing some hashing tasks with 6 cores, k-mers of length 33, and a hash size of 2^22:

	$ bash HashCounting.sh 6 33 22

Next, calculate the SVD and cluster the hased k-mers with 6 cores, a hash size of 2^22 and a cluster threshold of 0.8:

	$ bash KmerSVDClustering.sh 6 22 .8

Finally, partition the original reads (4 cores):

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

LSA has been written to be highly efficient in analyzing very large collections. For data sets larger than a few Gb, significant gains in wall time can be made by running in a cluster environment. In these cases, the process is essentially the same as what is outlined above. Detailed steps and job array submission scripts can be found in LSFscripts/README.md.