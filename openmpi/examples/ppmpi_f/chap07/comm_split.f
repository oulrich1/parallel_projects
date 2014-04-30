C  comm_split.f -- build a collection of q communicators using 
C      MPI_Comm_split
C 
C  Input: none
C  Output:  Results of doing a broadcast across each of the q 
C      communicators.
C 
C  Note:  Assumes the number of processes, p = q^2
C 
C  See Chap. 7, pp. 120 & ff in PPMPI
C 
      PROGRAM Split
      INCLUDE 'mpif.h'
      integer       p
      real          preal
      integer       my_rank
      integer       my_row_comm
      integer       my_row
      integer       q
      integer       test
      integer       my_rank_in_row
      integer       ierr
C
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      preal = p
      q =  sqrt(preal)
C
C  my_rank is rank in MPI_COMM_WORLD.
C  q*q = p   
      my_row = my_rank/q
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, my_row, my_rank,
     +      my_row_comm, ierr)
C
C  Test the new communicators   
      call MPI_COMM_RANK(my_row_comm,  my_rank_in_row, ierr)
      if (my_rank_in_row .EQ. 0) then
          test = my_row
      else
          test = 0
      endif
C
      call MPI_BCAST( test, 1, MPI_INTEGER, 0, my_row_comm, ierr)
C
      print *, 'Proc' ,my_rank, '> row= ',my_row,
     +      ', rank_in_row=',my_rank_in_row, ',test =', test
C
      call MPI_FINALIZE(ierr)
      end
