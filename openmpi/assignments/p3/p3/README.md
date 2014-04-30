just run with the command 
(the following contains -O2 Optimizations enabled)

mpicc -g -Wall -std=c99 \
                    -O2 -funroll-loops -march=native \
                    -o matrixmult ./matrixmult.c  \
                    -Wno-unused-result  \
                    -Wno-unused-variable


its also in the program code file

or just use ./make.sh (which just runs the above command)


then run with:

    mpirun -n 4 matrixmult < t1000.in > out



This program was run on jaguar and should compile
If it does not, please let me know if that is appropriate


If you would like to see the result of the matrix multiplication
change the flag "STDIN_RESULT_OUTPUT" from false to true

