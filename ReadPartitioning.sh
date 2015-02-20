###
# Dependencies
# numpy
# scipy
#

if [ "$#" -ne 1 ]; then
	echo "Illegal number of parameters"
	echo "Usage: bash KmerSVDClustering.sh numThreads"
	exit 1
fi

numThreads=$1

# ReadPartitions
mkdir tmp
numInputFiles=$(ls -l hashed_reads/*.hashq.* | grep ^- | wc -l)
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) partitioning reads in hashed input file {}; \
python LSA/write_partition_parts.py -r {} -i hashed_reads/ -o cluster_vectors/ -t tmp/ >> Logs/ReadPartitions.log 2>&1; \
if [ $? -ne 0 ]; then exit 1; fi' \
::: $(seq 1 $numInputFiles)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/ReadPartitions.log; exit 1; fi
rm -r tmp

# MergeIntermediatePartitions
mkdir read_partitions
numClusterTasks=`sed -n '1p' cluster_vectors/numClusters.txt`
parallel -j $numThreads --no-notice --halt-on-error 2 \
'echo $(date) merging partitions parts for cluster {}; \
python LSA/merge_partition_parts.py -r {} -i cluster_vectors/ -o read_partitions/ >> Logs/MergeIntermediatePartitions.log 2>&1; \
if [ $? -ne 0 ]; then exit 1; fi' \
::: $(seq 1 $numClusterTasks)
if [ $? -ne 0 ]; then echo "printing end of last log file..."; tail Logs/MergeIntermediatePartitions.log; exit 1; fi

echo $(date) Read partitioning is complete


echo Quantifying spike enrichment
spikeCount=`grep '@Spike' original_reads/*.fastq | wc -l`
echo Total spiked reads: $spikeCount
for i in $(seq 1 $numClusterTasks)
  do c=`grep '@Spike' read_partitions/$((i-1))/*.fastq | wc -l`
  t=`find read_partitions/$((i-1)) -type f -exec wc -l {} \; | awk '{total += $1} END{print total}'`
  echo partition $((i-1)): $c spiked reads out of $((t/4)) total reads in partition $((i-1)) \($(($c*400/$t))%\) >> read_partitions/spikeCoutns.txt
done

echo "Spiked read counts by partition (top 5)"
sort -k3,3nr read_partitions/spikeCoutns.txt | head -n5