#! /bin/sh

# 1, 2, 5, 10, 20, 30
# for each run in [0, 4]
#   execute all the following
# when finished, find min time

rundirectory="run5"
echo "$rundirectory/test1.out"
mkdir "$rundirectory"
qsub -o $rundirectory/test1.out -l select=1:ncpus=1   -l "walltime=2:00:00" ./test1.sh
qsub -o $rundirectory/test2.out -l select=1:ncpus=2   -l "walltime=2:00:00" ./test2.sh
qsub -o $rundirectory/test5.out -l select=1:ncpus=5   -l "walltime=2:00:00" ./test5.sh
qsub -o $rundirectory/test10.out -l select=1:ncpus=10 -l "walltime=2:00:00" ./test10.sh
qsub -o $rundirectory/test20.out -l select=1:ncpus=20 -l "walltime=2:00:00" ./test20.sh
qsub -o $rundirectory/test30.out -l select=1:ncpus=30 -l "walltime=2:00:00" ./test30.sh