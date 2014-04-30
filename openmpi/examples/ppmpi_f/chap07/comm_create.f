C  comm_create.f -- builds a communicator from the first q processes
C      in a communicator containing p = q^2 processes.
C 
C  Input: none
C  Output: q -- program tests correct creation of new communicator
C      by broadcasting the value 1 to its members -- all other 
C      processes have the value 0 -- global sum computed across
C      all the processes.
C
C  Note:  Assumes that MPI_COMM_WORLD contains p = q^2 processes
C 
C  See Chap 7, pp. 117 & ff in PPMPI
C 
      PROGRAM ComCrt
      INCLUDE 'mpif.h'
      integer        MAX_PROCS
      parameter      (MAX_PROCS = 100)
      integer        p
      real           p_real
      integer        q   
      integer        my_rank
      integer        group_world
      integer        first_row_group
      integer        first_row_comm
      integer        process_ranks(0:MAX_PROCS-1)
      integer        proc
      integer        test
      integer        sum
      integer        my_rank_in_first_row
      integer        ierr
C
C
      test = 0
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      p_real = p
      q =  sqrt(p_real)
C
C  Make a list of the processes in the new
C  communicator   
      do 100 proc = 0,  q-1  
          process_ranks(proc) = proc
 100  continue
C
C  Get the group underlying MPI_COMM_WORLD   
      call MPI_COMM_GROUP(MPI_COMM_WORLD,  group_world, ierr )
C
C  Create the new group   
      call MPI_GROUP_INCL(group_world, q, process_ranks,
     +                   first_row_group, ierr)
C
C  Create the new communicator   
      call MPI_COMM_CREATE(MPI_COMM_WORLD, first_row_group,
     +      first_row_comm, ierr)
C
C  Now check whether we can do collective ops in first_row_comm   
      if (my_rank .LT. q)  then
          call MPI_COMM_RANK(first_row_comm,  
     +                       my_rank_in_first_row, ierr)
          if (my_rank_in_first_row .EQ. 0) then
               test = 1
          endif
          call MPI_BCAST( test, 1, MPI_INTEGER, 0, 
     +                    first_row_comm, ierr)
      endif
      call MPI_REDUCE( test,  sum, 1, MPI_INTEGER, MPI_SUM, 0, 
     +                 MPI_COMM_WORLD, ierr )
C
      if (my_rank .EQ. 0)  then
          print *,'q = ', q, ' sum = ', sum
      endif
C
      call MPI_FINALIZE(ierr)
      end
