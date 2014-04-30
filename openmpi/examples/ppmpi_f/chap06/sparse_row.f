C  sparse_row.f -- pack a row of a sparse matrix and send from process 0
C      to process 1.  
C 
C  Input: none
C  Output: the row received by process 1.
C 
C  Notes:  
C      1. This program should only be run with 2 processes.  
C      2. Only the row of the matrix is created on both processes.
C 
C  See Chap. 6, pp. 104 & ff in PPMPI
C 
      PROGRAM SpaRow
      INCLUDE 'mpif.h'
      integer  HUGE
      parameter (HUGE = 100)
      integer   p
      integer   my_rank
      real      entries(10)
      integer   column_subscripts(10)
      integer   nonzeroes
      integer   position
      integer   row_number
      character buffer *100    
      integer   status(MPI_STATUS_SIZE)
      integer   ierr 
      integer   i
      data nonzeroes /10/
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      if (my_rank .EQ. 0)  then
C  Get the number of nonzeros in the row.   
C  Initialize entries and column_subscripts   
          do 100 i = 1, nonzeroes   
              entries(i)           = 2*i
              column_subscripts(i) = 3*i
 100      continue
C
C  Now pack the data and send   
          position = 1
          call MPI_PACK( nonzeroes, 1, MPI_INTEGER, buffer, HUGE,
     +         position, MPI_COMM_WORLD, ierr )
          call MPI_PACK( row_number, 1, MPI_INTEGER, buffer, HUGE,
     +          position, MPI_COMM_WORLD, ierr )
          call MPI_PACK(entries, nonzeroes, MPI_REAL , buffer,
     +         HUGE,  position, MPI_COMM_WORLD, ierr )
          call MPI_PACK(column_subscripts,nonzeroes,MPI_INTEGER,
     +         buffer, HUGE,  position, MPI_COMM_WORLD, ierr )
          call MPI_SEND(buffer, position, MPI_PACKED, 1, 0,
     +         MPI_COMM_WORLD, ierr )
       else    
          call MPI_RECV(buffer, HUGE, MPI_PACKED, 0, 0,
     +         MPI_COMM_WORLD,  status, ierr )
          position = 1
          call MPI_UNPACK(buffer, HUGE,  position,  nonzeroes,
     +         1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
          call MPI_UNPACK(buffer, HUGE,  position,  row_number,
     +         1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
          call MPI_UNPACK(buffer,HUGE,  position, entries,
     +         nonzeroes, MPI_REAL , MPI_COMM_WORLD, ierr )
          call MPI_UNPACK(buffer, HUGE,  position, 
     +         column_subscripts,
     +         nonzeroes, MPI_INTEGER, MPI_COMM_WORLD, ierr )
          do 200 i = 1, nonzeroes  
              print *, entries(i), column_subscripts(i)
 200      continue
      endif
C
      call MPI_FINALIZE(ierr)
      end
