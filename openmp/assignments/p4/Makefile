
CC = gcc 
CFLAGS= -g -Wall -fopenmp -std=c99 -fopenmp \
		-Wno-unused-result -Wno-unused-variable  \
		# -O2 -funroll-loops -march=native 
				# -pg -fprofile-arcs -ftest-coverage \
		 

LIBFLAGS= -lm

INCDIR= -I./

all: gaussian

gaussian: gaussian.o matrix.o 
	$(CC) $(CFLAGS)  $(INCDIR) $+ -o $@ $(LIBFLAGS)

gaussian.o: gaussian.c
	$(CC) $(CFLAGS)  $(INCDIR) -c $< $(LIBFLAGS)

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS)  $(INCDIR) -c $< $(LIBFLAGS)


clean:
	rm -rf *.o gaussian
