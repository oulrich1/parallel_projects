C  comm_test.f -- creates a communicator from the first q processes
C      in a communicator containing p = q^2 processes.  Broadcasts
C      an array to the members of the newly created communicator.
C 
C  Input: none
C  Output: Contents of array broadcast to each process in the newly
C      created communicator
C
C  Note:  MPI_COMM_WORLD should contain p = q^2 processes.
C 
C  See Chap 7., pp. 117 & ff. in PPMPI
C 
      PROGRAM Test
      INCLUDE 'mpif.h'
      integer        MAX_PROCS, n_bar
      parameter      (MAX_PROCS = 100, n_bar = 2)
      integer        p
      real           p_real
      integer        q   
      integer        my_rank
      integer        group_world
      integer        first_row_group
      integer        first_row_comm
      integer        process_ranks(0:MAX_PROCS-1)
      integer        proc, i
      real           A_00(n_bar * n_bar)
      integer        my_rank_in_first_row
      integer        ierr
C
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      p_real = p
      q = sqrt(p_real)
C
C  Make a list of the processes in the new
C  communicator   
      do 100 proc = 0, q-1
          process_ranks(proc) = proc
 100  continue
C
C  Get the group underlying MPI_COMM_WORLD   
      call MPI_COMM_GROUP(MPI_COMM_WORLD,group_world, ierr)
C
C  Create the new group   
      call MPI_GROUP_INCL(group_world, q, process_ranks,
     +      first_row_group, ierr)
C
C  Create the new communicator   
      call MPI_COMM_CREATE(MPI_COMM_WORLD, first_row_group,
     +      first_row_comm, ierr)
C
C  Now broadcast across the first row   
      if (my_rank .LT. q)  then
          call MPI_COMM_RANK(first_row_comm,  
     +                       my_rank_in_first_row, ierr)
C
          if (my_rank_in_first_row .EQ. 0)  then
C  Initialize A_00   
                do 200 i = 1 , n_bar*n_bar  
                     A_00(i) = i
 200            continue
          endif
          call MPI_BCAST(A_00, n_bar*n_bar, MPI_REAL , 0,
     +         first_row_comm, ierr)
C
          print *,'Process ', my_rank, ' > '
              print *,  A_00
          print *, ' '
      endif
      call MPI_FINALIZE(ierr)
      end
