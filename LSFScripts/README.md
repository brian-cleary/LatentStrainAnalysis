LSA has been optimized for the analysis of very large data sets in a highly distributed computing environment. To run this code you'll need Python (2.7), NumPy, SciPy, Gensim, and Pyro4. The following steps, which follow the same procedure as the streamlined "Getting Started" version but are more verbose, are generally one of two types: a command to create an LSF submission script, or a command to submit an LSF job.

To get started, download the LSA repo (eg via git clone). From here on, we'll assume that you've got the repo and your project housed in one location. You'll need to change the queue ("-q") argument in each *.q file to suit your environment.

Initialize your project with a few directories ("-i" takes the location of your reads, and "-n" takes the number of samples):

	$ python LSFScripts/setupDirs.py -i /input/reads/ -n 50
		
Begin by splitting the original reads (from many samples) into many small files

	In array_merge.py:
		
			This assumes the files are named sample_id.*.fastq.1 and sample_id.*.fastq.2 for paired reads. If you used some other naming convention, this needs to be reflected in line 26
			
	$ bsub < LSFScripts/SplitInput_ArrayJob.q
	
	### WARNING ###
	
		If your input files are significantly different from paired fastq files separated into 2 parts (.fastq.1 and .fastq.2) plus a singleton file (.single.fastq.1, then you will either need to modify these python files, or just take it upon yourself to split your files into chunks containing ~1million reads each, and named like: sample_id.fastq.xxx, where ".xxx" is the chunk number (eg '.021')
		

Create a k-mer hash function by drawing a bunch of random hyperplanes. If you want to adjust the k-mer length or hash size, alter the "-k" or "-s" arguments in the create_hash.py command of CreateHash_Job.q.

	$ python LSFScripts/create_jobs.py -j CreateHash -i ./
		
	$ bsub < LSFScripts/CreateHash_Job.q
		
	Look at the log files for this when it's done. Also, the hash function should be stored in hashed_reads/Wheels.txt
	

Hashing all the reads

	$ python LSFScripts/create_jobs.py -j HashReads -i ./
		
	$ bsub < LSFScripts/HashReads_ArrayJob.q
		
	Failure of a small fraction of these jobs is tolerable.
	

Tabulating k-mer counts in 1/5th of each sample

	First, make sure that there is just one *.fastq file per sample in original_reads/.
	
	The reason this is important is that the number of *.fastq files will be used to determine the array size. (The *.fastq.* files are no longer needed, so you can remove those as well if you want).
	
	$ python LSFScripts/create_jobs.py -j MergeHash -i ./
		
	$ bsub < LSFScripts/MergeHash_ArrayJob.q
		
	You don't really want any of these tasks to fail. So take a look at the logs when it's done and resubmit anything that died.
	

Merging the 5 count files for each sample

	$ python LSFScripts/create_jobs.py -j CombineFractions -i ./
		
	$ bsub < LSFScripts/CombineFractions_ArrayJob.q
		
	If any of these fail, they need to be run again.
	

Global (k-mer) conditioning

	$ python LSFScripts/create_jobs.py -j GlobalWeights -i ./
		
	$ bsub < LSFScripts/GlobalWeights_Job.q
		
	This launches a single job that must succeed to continue. Should produce cluster_vectors/global_weights.npy
	

Writing martix rows to separate files and local (sample) conditioning

	$ python LSFScripts/create_jobs.py -j KmerCorpus -i ./
		
	$ bsub < LSFScripts/KmerCorpus_ArrayJob.q
		
	These must all complete to continue. Relaunch any that failed. This job should produce one hashed_reads/*.conditioned file per sample.
	

Calculating the SVD (streaming!)

	$ python LSFScripts/create_jobs.py -j KmerLSI -i ./
	
	$ bsub < LSFScripts/KmerLSI_Job.q
		
	For very large matrices, this one will probably take a couple days to complete. Will produce cluster_vectors/kmer_lsi.gensim.
	

Create the cluster index

	$ python LSFScripts/create_jobs.py -j KmerClusterIndex -i ./
		
	$ bsub < LSFScripts/KmerClusterIndex_Job.q
		
	This step will set the k-mer cluster seeds, and the number of these seeds ultimately affects the resolution of partitioning. It is highly recommended that you check cluster_vectors/numClusters.txt for the number of clusters. If the resolution is markedly different from the expected / desired resolution, this job should be re-run with a different "-t" value in the submission script. Roughly speaking, we've found the following values to work for different scale datasets: 0.5-0.65 for large scale (Tb), 0.6-0.8 for medium scale (100Gb), >0.75 for small scale (10Gb). See misc/parameters.xlsx for more info.
	

Cluster blocks of k-mers

	$ bsub < LSFScripts/KmerClusterParts_ArrayJob.q
		

Merge cluster blocks

	$ bsub < LSFScripts/KmerClusterMerge_ArrayJob.q
		

Arrange k-mer clusters on disk

	$ python LSFScripts/create_jobs.py -j KmerClusterCols -i ./
		
	$ bsub < LSFScripts/KmerClusterCols_Job.q
		
	This should produce (among other things) a file cluster_vectors/kmer_cluster_sizes.npy
	

Partition all the read chunks

	$ python LSFScripts/create_jobs.py -j ReadPartitions -i ./
		
	You'll need to modify ReadPartitions_ArrayJob.q to contain your tmp directory of choice.
	
	$ sed 's/TMPDIR/\/your\/tmp\/dir/g' < LSFScripts/ReadPartitions_ArrayJob.q | bsub
		
	If a few of these fail, it's not super critical, but if a large number fail you'll want to resubmit them.
	

Merge the partition chunks

	$ python LSFScripts/create_jobs.py -j MergeIntermediatePartitions -i ./
		
	$ bsub < LSFScripts/MergeIntermediatePartitions_ArrayJob.q
		
	If any of these jobs fail you'll need to resubmit them.
	

If you've made it this far...good job! Your reads are now partitioned. Have at em'!