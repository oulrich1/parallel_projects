# Makefile for building program that solves a random dense linear system
#     using ScaLAPACK on a network of SGI's running mpich.  Change macro 
#     definitions to suit your system.

MPI_DIR = /usr/local/mpich
MPI_LIB_DIR = $(MPI_DIR)/lib/IRIX/ch_p4
SCALAPACK_DIR = /home/peter/SCALAPACK
BLACS_DIR = /home/peter/BLACS/LIB

F77 = f77
INCLUDE = 
# The "-u" tells sgi fortran to have no implicit typing
# The "-w0" tells sgi fortran to warn of unused variables
# The "-col72" tells sgi fortran to ignore code after column 72
F77_FLAGS = -u -g -w0 -col72

MPI_LIB = -L$(MPI_LIB_DIR) -lmpi
BLACS_LIB = $(BLACS_DIR)/blacsF77init_MPI-SGI-1.a \
	$(BLACS_DIR)/blacsCinit_MPI-SGI-1.a \
	$(BLACS_DIR)/blacs_MPI-SGI-1.a
LIB_BLAS = -lblas
LIB_PBLAS = $(SCALAPACK_DIR)/pblas_SGI.a
LIB_SCALAPACK = $(SCALAPACK_DIR)/scalapack_SGI.a \
	$(SCALAPACK_DIR)/tools_SGI.a

LIB = $(LIB_SCALAPACK) $(LIB_SCALAPACK) $(LIB_SCALAPACK) \
	$(LIB_PBLAS) $(LIB_PBLAS) \
	$(LIB_BLAS) $(LIB_BLAS) \
	$(BLACS_LIB) $(BLACS_LIB) $(MPI_LIB) -lm

linsolve: mpif.h linsolve.o
	$(F77) -o linsolve linsolve.o $(LIB)

mpif.h:
	rm -f mpif.h
	ln -s $(MPI_DIR)/include/mpif.h mpif.h

clean:
	rm -f linsolve.o linsolve

.f.o:
	$(F77) -c $(INCLUDE) $(F77_FLAGS) $*.f
