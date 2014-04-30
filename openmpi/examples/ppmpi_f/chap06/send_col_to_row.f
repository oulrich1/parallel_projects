C  send_col_to_row.f -- send column 1 of a matrix on process 0 to row 1
C      on process 1.
C 
C  Input: none
C  Output: The row received by process 1.
C 
C  Notes:  
C      1.  This program should only be run with 2 processes 
C      2.  Since fortran matrices are stored in column-major order,
C          this program uses the derived type on the receiving 
C          process rather than on the sending process
C 
C  See Chap 6., pp. 98 & ff in PPMPI
C 
      PROGRAM SndCR
      INCLUDE 'mpif.h'
      integer p
      integer my_rank
      real    A(10 ,10)
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      integer   row_mpi_t
      integer   i, j
C
      call MPI_INIT(ierr)
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
 100      continue
          call MPI_SEND( A(1,1), 10, MPI_REAL, 1, 0,
     +              MPI_COMM_WORLD, ierr )
      else 
          do 300 i = 1, 10
              do 400 j = 1, 10
                  A(i,j) =  0.0
 400          continue
 300      continue
          call MPI_RECV( A(1,1), 1, row_mpi_t , 0, 0,
     +              MPI_COMM_WORLD,  status, ierr )
          do 500 j = 1, 10
              print *, A(1,j)
 500      continue
      endif
C
      call MPI_FINALIZE(ierr)
      end
