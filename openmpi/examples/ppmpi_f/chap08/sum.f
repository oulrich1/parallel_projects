C  sum.f -- add two vectors using cyclic distribution of arrays.  Program
C      to illustrate use of cyclic_io functions.
C 
C  Input:
C      n:  order of vectors
C      x, y:  the vectors being added
C 
C  Output:
C      z: the sum vector
C 
C  Note:  Link with cio_mod.o and cyclic_io.o
C 
C  See Chap 8, pp. 170 & ff in PPMPI
C
      PROGRAM Sum
C  Header file for the basic I/O functions   
      INCLUDE 'ciof.h'
C  Header file for the cyclic array I/O functions   
      INCLUDE 'cyclic_iof.h'
      INCLUDE 'mpif.h'
      integer      x(0 : 8 + 1024 + 1024 -1)
      integer      y(0 : 8 + 1024 + 1024 -1)
      integer      z(0 : 8 + 1024 + 1024 -1)
      integer      n
      integer      io_comm
      integer      i
      integer      retval
      integer      ierr
C
      call MPI_INIT( ierr)
C
C  Build communicator for I/O   
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval = Cache_io_rank(MPI_COMM_WORLD, io_comm )
      if (retval   .EQ.   NO_IO_ATTR) then
          call MPI_ABORT(MPI_COMM_WORLD, -1, ierr )
      endif
C
C  Get n   
      retval = Cread(io_comm, 1, 'Enter array order  ',  n)
C
C  Initialize scalar members.  Calls    
C  function for building derived type   
      call INITialize_params( io_comm, n,  x)
      call INITialize_params( io_comm, n,  y)
      call INITialize_params( io_comm, n,  z)
C
C  Get vector elements   
      call Read_entries('Enter elements of x ',  x)
      call Read_entries('Enter elements of y ',  y)
C
C  Add local entries
      do 100 i = 0 , z(local_size)-1   
          z(local_entry + i) =
     +        x(local_entry + i) + y(local_entry +i)
      print *, ' sum (',i,') = ', z(local_entry+i)
      print *, x(local_entry + i), y(local_entry + i)
 100  continue
C
C  Print z   
      call Print_params(z)
      call Print_entries('x + y =',  z)
C
      call MPI_FINALIZE(ierr)
      end
