#!/bin/bash
#BSUB -J SplitInput[1-176]
#BSUB -o /project/home/Logs/SplitInput-Out-%I.out
#BSUB -e /project/home/Logs/SplitInput-Err-%I.err
#BSUB -q week
#BSUB -W 23:58
source /broad/software/scripts/useuse
reuse Python-2.7
export PYTHONPATH=/home/unix/bcleary/src/lib/python2.7/site-packages:$PYTHONPATH
cd /project/home
echo Date: `date`
t1=`date +%s`
sleep ${LSB_JOBINDEX}
python array_merge.py -r ${LSB_JOBINDEX} -i /input/data/ -o /project/home/original_reads/
[ $? -eq 0 ] || echo 'JOB FAILURE: $?'
echo Date: `date`
t2=`date +%s`
tdiff=`echo 'scale=3;('$t2'-'$t1')/3600' | bc`
echo 'Total time:  '$tdiff' hours'