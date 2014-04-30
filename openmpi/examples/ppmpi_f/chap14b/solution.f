C  solution.f -- functions for keeping track of the best solution found
C      so far -- for use in parallel tree search.
C 
C  See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
C      search.
C
C *******************************************************************  
      integer  function Best_solution_funct(comm )
      integer   comm   
C   
      include 'mainf.h'
      include 'solutionf.h'
      include 'mpif.h'
      integer   recv_status(MPI_STATUS_SIZE)
      integer   probe_status(MPI_STATUS_SIZE)
      integer   done
      logical   message_pending
      integer   temp , i
C
      done = FALSE
      do while (done .EQ.  FALSE) 
          call MPI_IPROBE(MPI_ANY_SOURCE, SOLUTION_TAG, comm,
     +          message_pending,  probe_status, ierr)
          if (message_pending)  then
              call MPI_RECV(temp_solution, solution_size, 
     +            MPI_INTEGER,
     +            probe_status(MPI_SOURCE), SOLUTION_TAG,
     +            comm, recv_status, ierr)
              if (temp_solution(0) .LT. best_solution(0))  then
                  do 100 i = 0, solution_size-1
                      temp  = temp_solution(i)
                      temp_solution(i) = best_solution(i)
                      best_solution(i) = temp
 100              continue
              endif
          else  
              done = TRUE
          endif
      end do
C
      Best_solution_funct = Local_best_solution( )
      return
      end   
C
C
C ********************************************************>**********  
      integer function Local_best_solution( ) 
C
      include 'solutionf.h'
      Local_best_solution = best_solution(0)
      return
      end
C
C ******************************************************************   
      subroutine Local_solution_update(cost, node)
      integer  cost   
      integer  node(0:15)
C
      include 'solutionf.h'
      include 'mainf.h'
      include 'node_stackf.h'
C   
      integer  i
      integer  bs_ptr
      integer  a_ptr
C
      best_solution(0) = cost 
      bs_ptr = 1
      a_ptr =  Ancestors
      do 100  i = 0,  max_depth -1
           best_solution(bs_ptr+i) = node(a_ptr+i)
100   continue
      best_solution(bs_ptr + i) = node(Sibling_rank)
      if (TREE_DEBUG) then
         print *, 'Local Update Soln =', bs_ptr,
     +             best_solution(bs_ptr)
      endif
C
      return
      end   
C
C
C ******************************************************************   
C  Assumes plenty of system buffering.  So beware!   
      subroutine Bcast_solution(comm)
      integer comm
C
      include 'mpif.h'
      include 'solutionf.h'
      include 'mainf.h'
      integer q
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
      call MPI_COMM_SIZE(comm,  p, ierr)
      do 100 q = 0 ,my_rank-1
      if (TREE_DEBUG) then
          print *, my_rank,' Send Best Soln'
      endif
          call MPI_SEND(best_solution, solution_size, 
     +        MPI_INTEGER, q,
     +        SOLUTION_TAG, comm, ierr)
 100  continue
      do 200 q = my_rank+1 , p -1
      if (TREE_DEBUG) then
          print *,my_rank, ' Send Best Soln'
      endif
          call MPI_SEND(best_solution, solution_size,
     +        MPI_INTEGER, q,
     +        SOLUTION_TAG, comm, ierr)
 200  continue
      return
      end  
C
C *******************************************************************  
C  Return 0 if OK, -1 otherwise   
      integer function Initialize_soln( )
C
      include 'mainf.h'
      include 'solutionf.h'
C
      integer i
C
      solution_size = max_depth + 2
C  
      best_solution(0) = INFINITY
C
      best_solution(1) = 0          
      do 200 i = 1, solution_size-1
          best_solution(i) = -1
          temp_solution(i) = -1
 200  continue
C
      Initialize_soln = 0
      return
      end
C
C *******************************************************************
C   Assumes each process has access to stdout
      subroutine Print_local_soln(comm)
      integer comm
C
      include 'mpif.h'
      include 'mainf.h'
      include 'solutionf.h'   
      integer i
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      print *, '*************************************************'
      print *, 'Process ',my_rank, ' > Minimum cost = ', 
     +     Local_best_solution()
      print *,'Min cost path = '
          print *, '   ',( best_solution(i), i = 1, solution_size-1)
      return
      end
C
C *******************************************************************  
      subroutine Print_solution(io_comm)
      integer io_comm   
      integer io_rank
      integer my_rank
      integer i
      integer retval, ierr
C
      include 'mpif.h'
      include 'solutionf.h'   
      include 'ciof.h'
C
      call MPI_COMM_RANK(io_comm,  my_rank, ierr)
      retval = Get_io_rank(io_comm,  io_rank)
      if (my_rank .EQ. io_rank)  then
          print *,'Minimum cost = ', Local_best_solution()
          print *, 'Path to minimum cost = '
          print *,'    ', (best_solution(i), i = 1, solution_size-1)
      endif
C
      return
      end   
C
C *******************************************************************  
C  Find owner of global best solution and broadcast   
      subroutine Update_solution(comm)
      integer  comm
      integer  local_data(0:1)
      integer  global_data(0:1)
      include 'mainf.h'
      include 'mpif.h'
      include 'solutionf.h'
C
      local_data(soln) = Local_best_solution()
      call MPI_COMM_RANK(comm,  local_data(rank),ierr)
      call MPI_ALLREDUCE( local_data,  global_data, 1, 
     +     MPI_2INTEGER,
     +     MPI_MINLOC, comm, ierr)
C
      if (TREE_DEBUG) then
         print *, ' Broadcast Best solution!, size = ',
     +   solution_size
      endif
      call MPI_BCAST(best_solution, solution_size, MPI_INTEGER,
     +     global_data(rank), comm, ierr)
      return
      end   
C
