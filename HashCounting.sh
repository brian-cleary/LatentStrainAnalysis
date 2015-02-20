###
# Dependencies
# numpy
# scipy
# gensim
# GNU parallel (citation info at bottom)
###

if [ "$#" -ne 3 ]; then
	echo "Illegal number of parameters"
	echo "Usage: bash HashCounting.sh numThreads kmerLen hashSize"
	exit 1
fi

numThreads=$1
kmerLen=$2
export hashSize=$3

mkdir Logs

# CreateHash
echo $(date) Creating the hash function with k-mers of length $kmerLen and hash size 2^$hashSize
mkdir hashed_reads
python LSA/create_hash.py -i original_reads/ -o hashed_reads/ -k 33 -s $hashSize > Logs/CreateHash.log 2>&1
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/CreateHash.log; exit 1; fi

# HashReads
numInputFiles=$(ls -l original_reads/*.fastq | grep ^- | wc -l)
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) hashing reads in file {}; \
python LSA/hash_fastq_reads.py -r {} -i original_reads/ -o hashed_reads/ >> Logs/HashReads.log 2>&1' \
::: $(seq 1 $numInputFiles)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/HashReads.log; exit 1; fi

# MergeHash / CombineFractions
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) counting hashed k-mers in file {}; \
python LSA/merge_hashq_files.py -r {} -i hashed_reads/ -o hashed_reads/ >> Logs/MergeHash.log 2>&1; \
python LSA/merge_hashq_fractions.py -r {} -i hashed_reads/ -o hashed_reads/ >> Logs/CombineFractions.log 2>&1' \
::: $(seq 1 $numInputFiles)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/MergeHash.log; tail Logs/CombineFractions.log; exit 1; fi

# GlobalWeights
echo $(date) Finding global weights for each hashed k-mer
mkdir cluster_vectors
python LSA/tfidf_corpus.py -i hashed_reads/ -o cluster_vectors/ > Logs/GlobalWeights.log 2>&1
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/GlobalWeights.log; exit 1; fi

# KmerCorpus
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) writing k-mer corpus for file {}; \
python LSA/kmer_corpus.py -r {} -i hashed_reads/ -o cluster_vectors/ >> Logs/KmerCorpus.log 2>&1' \
::: $(seq 1 $numInputFiles)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/KmerCorpus.log; exit 1; fi

echo $(date) HashCounting is finished



# {Tange2011a,
# title = {GNU Parallel - The Command-Line Power Tool},
# author = {O. Tange},
# address = {Frederiksberg, Denmark},
# journal = {;login: The USENIX Magazine},
# month = {Feb},
# number = {1},
# volume = {36},
# url = {http://www.gnu.org/s/parallel},
# year = {2011},
# pages = {42-47}
# }