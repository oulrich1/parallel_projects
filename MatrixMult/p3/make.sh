#! /usr/bin/env bash
echo "mpicc -g -Wall -std=c99 \
                    -funroll-loops -march=native \
                    -o matrixmult ./matrixmult.c  \
                    -Wno-unused-result  \
                    -Wno-unused-variable -lm "
mpicc -g -Wall -std=c99 \
                    -funroll-loops -march=native \
                    -o matrixmult ./matrixmult.c  \
                    -Wno-unused-result  \
                    -Wno-unused-variable -lm 

