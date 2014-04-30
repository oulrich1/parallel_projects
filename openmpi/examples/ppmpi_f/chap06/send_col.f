C  send_col.f -- send the third column of a matrix from process 0 to
C      process 1
C 
C  Input:  None
C  Output:  The column received by process 1
C 
C  Notes:  
C      1.  This program should only be run with 2 processes
C      2.  Since matrices in fortran are stored in column-major order,
C          we can send a column without using a derived datatype.
C          Program send_row.f uses MPI_Type_vector to send a row.
C 
C  See Chap 6., p. 96, in PPMPI
C 
      PROGRAM SenCol
      INCLUDE 'mpif.h'
      integer   p
      integer   my_rank
      real      A(10 ,10)
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      integer i, j
C
      call MPI_INIT( ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      if (my_rank .EQ. 0)  then
          do 100 i = 1, 10   
              do 200 j =1, 10 
                  A(i, j) = j
 200          continue
 100      continue
          call MPI_SEND( A(1,3), 10, MPI_REAL , 1, 0,
     +           MPI_COMM_WORLD, ierr )
      else    
          call MPI_RECV( A(1,3 ), 10, MPI_REAL , 0, 0,
     +         MPI_COMM_WORLD,  status, ierr )
          do 300  i = 1, 10
              print*, A(i, 3 )
 300      continue
      endif
C
      call MPI_FINALIZE(ierr)
      end
