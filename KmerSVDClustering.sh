###
# Dependencies
# numpy
# scipy
# gensim: pip install gensim
# Pyro4: pip install Pyro4
###

if [ "$#" -ne 3 ]; then
	echo "Illegal number of parameters"
	echo "Usage: bash KmerSVDClustering.sh numThreads hashSize clustThresh"
	exit 1
fi

numThreads=$1
hashSize=$2
export clustThresh=$3

# KmerLSI
echo $(date) Starting streaming SVD of conditioned k-mer abundance matrix
python LSA/kmer_lsi.py -i hashed_reads/ -o cluster_vectors/ -s > Logs/KmerLSI.log 2>&1
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/KmerLSI.log; exit 1; fi

# KmerClusterIndex
echo $(date) Finding eigengenome seeds for clustering using a cosine distance threshold of $clustThresh
python LSA/kmer_cluster_index.py -i hashed_reads/ -o cluster_vectors/ -t $clustThresh > Logs/KmerClusterIndex.log 2>&1
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/KmerClusterIndex.log; exit 1; fi

# KmerClusterParts
echo $(date) Starting k-mer clustering
# number of tasks is 2**hash_size/10**6 + 1
numClusterTasks=`sed -n 1p hashed_reads/hashParts.txt`
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) clustering k-mer chunk {}; \
python LSA/kmer_cluster_part.py -r {} -i hashed_reads/ -o cluster_vectors/ -t $clustThresh >> Logs/KmerClusterParts.log 2>&1' \
::: $(seq 1 $numClusterTasks)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/KmerClusterParts.log; exit 1; fi

# KmerClusterMerge
# number of tasks is number of clusters
numClusterTasks=`sed -n '1p' cluster_vectors/numClusters.txt`
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) merging chunks for k-mer cluster {}; \
python LSA/kmer_cluster_merge.py -r {} -i cluster_vectors/ -o cluster_vectors/ >> Logs/KmerClusterMerge.log 2>&1' \
::: $(seq 1 $numClusterTasks)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/KmerClusterMerge.log; exit 1; fi

# KmerClusterCols
echo $(date) Arranging k-mer clusters on disk
python LSA/kmer_cluster_cols.py -i hashed_reads/ -o cluster_vectors/ > Logs/KmerClusterCols.log 2>&1
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/KmerClusterCols.log; exit 1; fi

echo $(date) Kmer clustering is complete