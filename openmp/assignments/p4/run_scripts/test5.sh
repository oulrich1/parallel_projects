#!/bin/sh 

##PBS -N myjob 
##PBS -j oe 
#export OMP_NUM_THREADS=10
$HOME/p4/gaussian 8000 5

