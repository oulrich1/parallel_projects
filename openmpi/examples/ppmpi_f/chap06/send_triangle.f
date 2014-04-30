C  send_triangle.f -- send the upper triangle of a matrix from process 0
C      to process 1
C 
C  Input:  None
C  Output:  The matrix received by process 1
C 
C  Note:  This program should only be run with 2 processes.
C 
C  See Chap 6, p. 98, in PPMPI.
C 
      PROGRAM SndTri
      INCLUDE 'mpif.h'
      integer   p
      integer   my_rank
      integer   A(10, 10)            
      integer   T(10, 10)            
      integer   displacements(0:9)
      integer   block_lengths(0:9)
      integer   index_mpi_t
      integer   i, j
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
C
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      do 100 i = 0,9 
          block_lengths(i) = 10-i
          displacements(i) = (10+1)*i
 100  continue
      call MPI_TYPE_INDEXED(10, block_lengths, displacements,
     +                         MPI_INTEGER ,  index_mpi_t, ierr)
      call MPI_TYPE_COMMIT( index_mpi_t, ierr)
C
      if (my_rank .EQ. 0)  then
          do 200 i = 1, 10  
              do 300 j = 1, 10
                  A(i ,j) =  i + j
 300          continue
 200      continue
          call MPI_SEND(A,1, index_mpi_t,1,0, MPI_COMM_WORLD, ierr)
       else  
          do 400 i = 1, 10
              do 500 j = 1, 10
                  T(i, j) = 0
 500          continue
 400      continue
          call MPI_RECV(T, 1, index_mpi_t, 0, 0, MPI_COMM_WORLD,  
     +                       status, ierr )
          do 600 i = 1, 10
            write (6,700) (T(i, j), j = 1,10)
 700        format(10I3)
 600      continue
      endif
C
      call MPI_FINALIZE(ierr)
      end
