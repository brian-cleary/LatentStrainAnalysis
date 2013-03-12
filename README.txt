Make a home for your project
$ mkdir /project/home
Put your input fastq files in the project
$ mkdir /project/home/original_reads
$ mv *.fastq /project/home/original_reads

Throughout the analysis, additional directories will be created:
hashed_reads/
cluster_vectors/
cluster_vectors/intermediate_clusters/	<--can be removed after merge_read_clusters.py
read_partitions/

READ PAIRS ARE EXPECTED TO BE INTERLEAVED
$ python merge_read_pairs.py trimmed.1.fastq trimmed.2.fastq pairs.fastq

SPLIT THE INPUT FILES TO ALLOW FOR GREATER DISTRIBUTION OF JOBS
$ python split_fastq_files.py -i /project/home/original_reads -s 1000

CREATE THE HASH
$ mkdir /project/home/hashed_reads
$ python create_hash.py -i /project/home/original_reads/ -o /project/home/hashed_reads/ -k kmer_length -s hash_size

HASH INPUT READS, CREATE HASH COUNTS
$ bsub < HashReads_ArrayJob.q
$ bsub < MergeHash_ArrayJob.q

CONDITION, SVD, AND CLUSTER KMER ABUNDANCE MATRIX
$ mkdir /project/home/cluster_vectors
$ bsub < KmerClusters_Job.q

WRITE VERBOSE READS INTO CLUSTERS
$ mkdir /project/home/cluster_vectors/intermediate_clusters
$ bsub < ReadPartitions_ArrayJob.q
$ mkdir /project/home/read_partitions
$ bsub < MergeIntermediatePartitions_Job.q
# remove intermediate_clusters/ if all looks good.

ANALYZE CLUSTER CONTENT
$ bsub < PartitionAnalysis_ArrayJob.q
# find jobs that died
$ egrep -ir "job killed" Logs/TongueDorsumAssemblyIDBA-Out-*