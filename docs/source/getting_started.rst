Getting Started
===============

Dependencies
^^^^^^^^^^^^

LSA depends on the following software and libraries::

	Python (2.7)
	NumPy
	SciPy
	Gensim
	GNU Parallel
	Pyro4

To install GNU Parallel (http://www.gnu.org/software/parallel/)::
	
	(wget -O - pi.dk/3 || curl pi.dk/3/) | bash
	
If you're having trouble with Parallel, check the version (parallel --version), and look for the following::

	WARNING: YOU ARE USING --tollef. IF THINGS ARE ACTING WEIRD USE --gnu.

In this case, heed the advice and use --gnu.

Note that you will not need Pyro4 to analyze the test data as described below.

Test Data
^^^^^^^^^

To get started, make sure you have the dependencies listed above (you'll need GNU Parallel, but not Pyro4), and then download the repo (eg via git clone).

Next, unpack testData.tar.gz (into original_reads/). This folder contains a subsample of 10k reads from each of the 18 metagenomic libraries used in Sharon et. al. (SRA052203), and in the original LSA methods paper.

Each sample has been randomly spiked with 0-2,000 mock reads from a *Bacillus thuringiensis* plasmid (NG_035027.1). Note that there is one file per sample, that mate pairs are interleaved, and that files are named sample_id.*.fastq. You'll need a couple GBs of RAM to work through the test data. If something goes wrong, it's best to first check the most recently modified file in "Logs/" to track down the error.

Begin by performing some hashing tasks using 6 threads, k-mers of length 33, and a hash size of 2^22:::

	$ bash HashCounting.sh 6 33 22

Next, calculate the SVD and cluster the hashed k-mers using 6 threads, a hash size of 2^22 and a cluster threshold of 0.8:::

	$ bash KmerSVDClustering.sh 6 22 .8

Finally, partition the original reads (4 threads):::

	$ bash ReadPartitioning.sh 4

Since the process includes several points of randomization, the results from each run will vary. However, when the read partitioning step completes, a report is made to indicate the specificity and completeness of spiked read partitioning:::

    Quantifying spike enrichment
    Total spiked reads: 21528
    Spiked read counts by partition (top 5)
    partition 5: 19771 spiked reads out of 20291 total reads in partition 5 (97%)
    partition 146: 1207 spiked reads out of 1557 total reads in partition 146 (77%)
    partition 131: 356 spiked reads out of 670 total reads in partition 131 (53%)
    partition 112: 142 spiked reads out of 492 total reads in partition 112 (28%)
    partition 129: 4 spiked reads out of 506 total reads in partition 129 (0%)

Analyzing larger collections
----------------------------

LSA has been written to be highly efficient in analyzing very large collections. For data sets larger than a few Gb, significant gains in wall time can be made by running in a cluster environment. In these cases, the process is essentially the same as what is outlined above. Detailed steps and job array submission scripts can be found in LSFScripts/README.md.
