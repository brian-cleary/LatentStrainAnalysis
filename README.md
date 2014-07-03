LSA was written to run in an LSF environment. Python 2.7 is required, and you'll need NumPy, SciPy, Gensim, and Pyro4 in your python environment. The steps below, which take you through LSA, are generally one of two types: a command to create an LSF submission script, or a command to submit an LSF job.

### WARNING ###
	Those who are not comfortable with LSF usage will likely find it very difficult to run this code. This code is not in the state of an out-of-the-box, command line utility. You will need to modify LSF submission scripts, check log files, and generally be able to tell when tasks have successfully completed.

Before getting started, place the code for LSA in some directory (henceforth "the code directory"), and create a home for your project (henceforth "/project/home"). You'll need to create the following subdirs in /project/home:

		$ mkdir /project/home/Logs
		
		$ mkdir /project/home/original_reads
		
		$ mkdir /project/home/hashed_reads
		
		$ mkdir /project/home/cluster_vectors
		
		$ mkdir /project/home/read_partitions
		

Begin by splitting the original reads (from many samples) into many small files

	Copy array_merge.py, SplitInput_ArrayJob.q, and merge_and_split_pair_files.py into /project/home
	
		In SplitInput_ArrayJob.q change:
		
			The size of the job array needs to be the number of samples in your collection (ie [1-50] for 50 samples)
			
			The paths to the -o and -e logs need to be your /project/home/Logs/
			
			Change line 10 to your project home: cd /project/home
			
			In line 14, the input ( -i ) needs to be the directory with the raw reads and the output ( -o ) needs to be /project/home/original_reads/
			
		In array_merge.py:
		
			This assumes the files are named sample_id.*.fastq.1 and sample_id.*.fastq.2 for paired reads. If you used some other naming convention, this needs to be reflected in line 26
			
	$ bsub < SplitInput_ArrayJob.q
	
	### WARNING ###
	
		If your input files are significantly different from paired fastq files separated into 2 parts (.fastq.1 and .fastq.2) plus a singleton file (.single.fastq.1, then you will either need to modify these python files, or just take it upon yourself to split your files into chunks containing ~1million reads each, and named like: sample_id.fastq.xxx, where ".xxx" is the chunk number (eg '.021')
		

Creating the hash function

	From the directory where you've placed LSA code:
	
		$ python create_jobs.py -j CreateHash -i /project/home/
		
	From /project/home:
	
		$ mkdir hashed_reads
		
		$ bsub < CreateHash_Job.q
		
	Look at the log files for this when it's done. Also, there should a file hashed_reads/Wheels.txt
	

Hashing all the reads

	From the code directory:
	
		$ python create_jobs.py -j HashReads -i /project/home/
		
	From /project/home
	
		$ bsub < HashReads_ArrayJob.q
		
	Failure of a small fraction of jobs is tolerable.
	

Tabulating k-mer counts in 1/5th of each sample

	First, make sure that there is just one *.fastq file per sample in original_reads/. Also, remove original_reads/random_kmers.fastq
	
	The reason this is important is that the number of *.fastq files will be used to determine the array size. (The *.fastq.* files are no longer needed, so you can remove those as well if you want).
	
	From the code directory:
	
		$ python create_jobs.py -j MergeHash -i /project/home/
		
	From /project/home:
	
		$ bsub < MergeHash_ArrayJob.q
		
	You don't really want any of these tasks to fail. So take a look at the logs when it's done and resubmit anything that died.
	

Merging the 5 count files for each sample

	From the code directory:
	
		$ python create_jobs.py -j CombineFractions -i /project/home/
		
	From /project/home:
	
		$ bsub < CombineFractions_ArrayJob.q
		
	If any of these fail, they need to be run again.
	

Global (k-mer) conditioning

	From the code directory:
	
		$ python create_jobs.py -j GlobalWeights -i /project/home/
		
	From /project/home:
	
		$ mkdir cluster_vectors
		
		$ bsub < GlobalWeights_Job.q
		
	This launches a single job that must succeed to continue. Should produce cluster_vectors/global_weights.npy
	

Writing martix rows to separate files and local (sample) conditioning

	From the code directory:
	
		$ python create_jobs.py -j KmerCorpus -i /project/home/
		
	From /project/home:
	
		$ bsub < KmerCorpus_ArrayJob.q
		
	These must all complete to continue. Relaunch any that failed. This job should produce one hashed_reads/*.conditioned file per sample.
	

Calculating the SVD

	From the code directory:
	
	
		$ python create_jobs.py -j KmerLSI_Job.q -i /project/home/
	From /project/home:
	
		$ bsub < KmerLSI_Job.q
		
	For very large matrices, this one will probably take a couple days to complete. Will produce cluster_vectors/kmer_lsi.gensim.
	

Create the cluster index

	From the code directory:
	
		$ python create_jobs.py -j KmerClusterIndex -i /project/home/
		
	From /project/home:
	
		$ bsub < KmerClusterIndex_Job.q
		
	This step will set the k-mer cluster seeds, and the number of these seeds ultimately affects the resolution of partitioning. It is highly recommended that you check the log file from this job for the number of clusters. If the resolution is markedly different from the expected / desired resolution, this job should be re-run with a different "-t" value in the submission script. Roughly speaking, we've found the following values to work for different scale datasets: 0.5-0.65 for large scale (Tb), 0.6-0.8 for medium scale (100Gb), >0.75 for small scale (10Gb).
	

Cluster blocks of k-mers

	From the code directory:
	
		$ python create_jobs.py -j KmerClusterParts -i /project/home/
		
	From /project/home:
	
	You will need to modify KmerClusterParts_ArrayJob.q so that it has the appropriate array size for the hash size we've chosen. The array size should be 2^hash_size/10^6 + 1, so if you stuck with the default of 31, this is 2^31/10^6 + 1 = 2148 -> [1-2148]
	
		$ bsub < KmerClusterParts_ArrayJob.q
		

Merge cluster blocks

	From the code directory:
	
		$ python create_jobs.py -j KmerClusterMerge -i /project/home/
		
	From /project/home/:
	
	You will need to modify KmerClusterMerge_ArrayJob.q so that it has the appropriate array size for the number of clusters. You can find the number of clusters printed at the end of the log KmerClusterIndex-Out.out. There will be a line like "cluster index has shape (2345,67), where the first of these numbers denotes the number of clusters.
	
		$ bsub < KmerClusterMerge_ArrayJob.q
		

Arrange k-mer clusters on disk

	From the code directory:
	
		$ python create_jobs.py -j KmerClusterCols -i /project/home/
		
	From /project/home/:
	
		$ bsub < KmerClusterCols_Job.q
		
	This should produce (among other things) a file cluster_vectors/kmer_cluster_sizes.npy
	

Partition all the read chunks

	From the code directory:
	
		$ python create_jobs.py -j ReadPartitions -i /project/home/
		
	From /project/home:
	
	You'll need to modify ReadPartitions_ArrayJob.q to contain your tmp directory of choice. Then replace the write_partition_parts.py option "-t TMPDIR" with "-t /your/tmp/dir".
	
		$ mkdir read_partitions
		
		$ bsub < ReadPartitions_ArrayJob.q
		
	If a few of these fail, it's not super critial, but if a large number fail you'll want to resubmit them.
	

Merge the partition chunks

	From the code directory:
	
		$ python create_jobs.py -j MergeIntermediatePartitions -i /project/home/
		
	From /project/home/:
	
		$ bsub < MergeIntermediatePartitions_ArrayJob.q
		
	If any of these jobs fail you'll need to resubmit them.
	

If you've made it this far...good job! Your reads are now partitioned. Have at em'.