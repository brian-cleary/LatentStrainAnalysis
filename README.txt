Make a home for your project
$ mkdir /project/home
Put your input fastq files in the project
$ mkdir /project/home/original_reads
$ mv *.fastq /project/home/original_reads

Throughout the analysis, additional directories will be created:
Logs/
hashed_reads/
cluster_vectors/
cluster_vectors/intermediate_clusters/	<--can be removed after merge_read_clusters.py
read_partitions/

READ PAIRS ARE EXPECTED TO BE INTERLEAVED
$ python merge_read_pairs.py trimmed.1.fastq trimmed.2.fastq pairs.fastq
INPUT FILES SHOULD HAVE NAMES LIKE:
SampleName.*.fastq
SampleName.trimmed.fastq
SampleName.lib1.trimmed.fastq
SampleName.lib2.trimmed.fastq
# data with the same "SampleName" prefix will be pooled for the purpose of finding the SVD, but will remain separate for assembling with variable library sizes.

SPLIT THE INPUT FILES TO ALLOW FOR GREATER DISTRIBUTION OF JOBS
$ python split_fastq_files.py -i /project/home/original_reads -s 1000

$ mkdir /project/home/Logs

# Default process:
SEMI-SUCCINCT PIPELINE
$ python doHash.py -i /project/home -k 41 -s 30
$ python doSVD.py -i /project/home
$ python doPartitions.py -i /project/home
$ python doMerge.py -i /project/home


# Step-by-step process:
VERBOSE PIPELINE

CREATE THE HASH
$ mkdir /project/home/hashed_reads
$ python create_hash.py -i /project/home/original_reads/ -o /project/home/hashed_reads/ -k kmer_length -s hash_size

HASH INPUT READS, CREATE HASH COUNTS
$ python create_jobs.py -j HashReads -i /project/home/
$ bsub < HashReads_ArrayJob.q
$ python create_jobs.py -j MergeHash -i /project/home/
$ bsub < MergeHash_ArrayJob.q

CONDITION, SVD, AND CLUSTER KMER ABUNDANCE MATRIX (streaming)
$ mkdir /project/home/cluster_vectors
$ python create_jobs.py -j GlobalWeights -i /project/home/
$ bsub < GlobalWeights_Job.q
$ python create_jobs.py -j LSIKmerClusters -i /project/home/
$ bsub < LSIKmerClusters_Job.q
CONDITION, SVD, AND CLUSTER KMER ABUNDANCE MATRIX (in memory)
$ mkdir /project/home/cluster_vectors
$ bsub < KmerClusters_Job.q

WRITE VERBOSE READS INTO CLUSTERS
$ python create_jobs.py -j ReadPartitions -i /project/home/
$ bsub < ReadPartitions_ArrayJob.q
$ mkdir /project/home/read_partitions
$ python create_jobs.py -j MergeIntermediatePartitions -i /project/home
$ bsub < MergeIntermediatePartitions_ArrayJob.q

ANALYZE CLUSTER CONTENT
$ bsub < PartitionAnalysis_ArrayJob.q
# find jobs that died
$ egrep -ir "job killed" Logs/TongueDorsumAssemblyIDBA-Out-*