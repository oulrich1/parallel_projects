
CC = gcc 
CFLAGS= -g -Wall -fopenmp -std=c99  		\
		# -O3 -funroll-loops -march=native 	\
		-Wno-unused-result -Wno-unused-variable 

LIBFLAGS= -lm

INCDIR= -I./

all: omp_hello

omp_hello: omp_hello.o omp_trap1.o omp_trap2b.o omp_trap3.o omp_workshare1.o
	$(CC) $(CFLAGS)  $(INCDIR) omp_hello.o 		-o omp_hello $(LIBFLAGS)
	$(CC) $(CFLAGS)  $(INCDIR) omp_trap1.o 		-o omp_trap1 $(LIBFLAGS)
	$(CC) $(CFLAGS)  $(INCDIR) omp_trap2b.o 	-o omp_trap2b $(LIBFLAGS)
	$(CC) $(CFLAGS)  $(INCDIR) omp_trap3.o 		-o omp_trap3 $(LIBFLAGS)
	$(CC) $(CFLAGS)  $(INCDIR) omp_workshare1.o  -o omp_workshare1 $(LIBFLAGS)

omp_hello.o: omp_hello.c
	$(CC) $(CFLAGS)  $(INCDIR)				-c omp_hello.c $(LIBFLAGS)

omp_trap1.o: omp_trap1.c
	$(CC) $(CFLAGS)  $(INCDIR)				-c omp_trap1.c $(LIBFLAGS)

omp_trap2b.o: omp_trap2b.c
	$(CC) $(CFLAGS)  $(INCDIR)				-c omp_trap2b.c $(LIBFLAGS)

omp_trap3.o: omp_trap3.c
	$(CC) $(CFLAGS)  $(INCDIR) 				-c omp_trap3.c $(LIBFLAGS)

omp_workshare1.o: omp_workshare1.c
	$(CC) $(CFLAGS)  $(INCDIR) 				-c $< $(LIBFLAGS)

clean:
	rm -rf *o 
	rm -rf omp_hello
	rm -rf omp_trap1
	rm -rf omp_trap2b
	rm -rf omp_trap3