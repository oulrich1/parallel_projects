#!/usr/bin/env bash
if [ -z $1 ]
then
    echo "Requires target name to build...exiting instead.."
    exit
fi
target=$1
gcc -c $target.c $(gsl-config --libs)
gcc $(gsl-config --libs) $target.o 

# openMP
# compile with -fopenmp
