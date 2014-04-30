C  send_row.f -- send third row of a matrix from process 0 to process 1
C  
C  Input: none
C  Output: the row received by process 1
C 
C  Notes:  
C      1.  This program should only be run with 2 processes
C      2.  Since matrices in fortran are stored in column-major order,
C          we use a derived datatype to send a row.  Program send_col.f 
C          sends a column without a derived type.
C 
C See Chap 6, pp. 96 & ff.,  in PPMPI
C
      PROGRAM SendRow
      INCLUDE 'mpif.h'
      integer   p
      integer   my_rank
      real      A(10,10)
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      integer   row_mpi_t
      integer i, j
C
      call MPI_INIT( ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      call MPI_TYPE_VECTOR(10, 1, 10, MPI_REAL, row_mpi_t, ierr)
      call MPI_TYPE_COMMIT( row_mpi_t, ierr)
C
      if (my_rank .EQ. 0)  then
          do 100 i = 1, 10  
              do 200 j = 1, 10 
                  A(i,j) =  i
 200          continue
 100       continue
           call MPI_SEND( A(3, 1), 1, row_mpi_t, 1, 0,
     +                  MPI_COMM_WORLD, ierr )
      else    
          call MPI_RECV( A(3, 1), 1, row_mpi_t, 0, 0,
     +         MPI_COMM_WORLD,  status, ierr )
          do 300 i = 1,10 
              print *, A(3, i)
 300      continue
      endif
C
      call MPI_FINALIZE(ierr)
      end
