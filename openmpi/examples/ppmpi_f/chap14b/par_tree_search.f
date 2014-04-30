C  par_tree_search.f -- controlling program for parallel tree search.
C      Process 0 generates initial tree and distributes it among the
C      processes.  Then each process runs through the basic loop
C 
C          do {
C              Local depth-first search for a while;
C              Handle requests for work;
C          } while (There's any work left on all processes);
C 
C  After exiting the loop, the global best solution is updated and
C  printed.
C  
C  See Chap 14, pp. 328 & ff., in PPMPI.
C
C
C *******************************************************************
      subroutine Par_tree_search(root, io_comm)
      integer       root      
      include  'mpif.h'
      include  'mainf.h'
      include  'statsf.h'
      include  'node_stackf.h'
      include  'work_remainsf.h'
C
      integer   node(0:15)
C
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr)
C
C  Generate initial set of nodes, 1 per process
      retval =  TRUE
      if (my_rank .EQ. 0)  then 
          call Generate(root, p,
     +                      MPI_COMM_WORLD, p, my_rank)
      endif
C
      call Scatter( node, MPI_COMM_WORLD)
C
      call Initialize(node, p, my_rank)
C
      retval = TRUE
      do while (retval .NE. FALSE)
C  Search for a while   
          if ( STATS) then
              call Start_time_rtn()
          endif
      if (TREE_DEBUG) then
         print *, my_rank, 'In PAR_TREE call PAR_DFS'
      endif
          call Par_dfs( MPI_COMM_WORLD, p, my_rank)
          if ( STATS ) then
               call Finish_time(par_dfs_time)
          endif
C
C  Service requests for work.   
          if ( STATS ) then
              call Start_time_rtn()
          endif
      if (TREE_DEBUG) then
           print *, my_rank, 'Call to Service Requests'
      endif
           call Service_requests(  MPI_COMM_WORLD)
      if (TREE_DEBUG) then
           print *, my_rank, 'Return from Service Requests'
      endif
          if ( STATS ) then
              call Finish_time(svc_req_time)
          endif
C
C  If local_stack isn't empty, return.            
C  If local_stack is empty, send
C  receive a message terminating program.
	      retval = Work_remains(MPI_COMM_WORLD, p, my_rank)
       if (TREE_DEBUG) then
          print *,my_rank, 'Return from Work_remains = ',
     +                   retval
       endif
       end do
C
C  Get global best solution
      if (TREE_DEBUG) then
          print *, my_rank, 'DONE WITH LOOP in PAR TREE
     +                       - Update Soln'
      endif
      call Update_solution(MPI_COMM_WORLD)
      call Print_solution(io_comm)
C
      return 
      end
C
C
C *******************************************************************  
C  Called in Par_tree_search   
      subroutine Generate(root,  size_in, comm,
     +                   p, my_rank)
      include 'par_tree_searchf.h'
      include 'par_dfsf.h'
      include 'mainf.h'
      include 'node_stackf.h'
      include 'statsf.h'
      include 'solutionf.h'
      include 'mpif.h'
C
      integer       root
      integer       size_in        
      integer      comm        
C
      integer       stack_size 
      integer       node(0:15)
      integer       temp_sol
C
      stack_size = 1
      do while ( Empty(local_stack) .NE. TRUE .AND.
     +          (stack_size .LT. size_in))
          call Pop(local_stack, node)
          if (STATS) then
              call Incr_stat(nodes_expanded)
          endif
C
          stack_size = stack_size - 1
          if (Solution(node) .NE. FALSE )  then
      if (TREE_DEBUG) then
          print *,my_rank, 'SOLUTION = TRUE'
      endif
              temp_sol = Local_evaluate(node)
              if (temp_sol .LT. Local_best_solution())  then
                  call Local_solution_update(temp_sol, node)
              endif
          else 
             if (Feasible(node, comm).NE. FALSE )  then
              call PTS_expand(node, stack_size, p, my_rank)
             endif 
          endif
      enddo
C
      if (stack_size .LT. size_in)  then
          print *, 'Too few nodes in tree!'
          print *, 'Aborting.'
          call MPI_ABORT(MPI_COMM_WORLD, -1, ierr )
      endif
C
      return
      end  
C
C
C *******************************************************************  
      integer function Local_evaluate(node)
      integer node(0:15)
      include 'node_stackf.h'
C
      Local_evaluate =  My_rand(node(Seed))
      return
      end   
C
C
C *******************************************************************
      subroutine PTS_expand(parent, stack_size, p, my_rank)
      integer     parent(0:15)               
      integer     stack_size
C
      include 'par_dfsf.h'
      include 'node_stackf.h'
      include 'mainf.h'
      include 'mpif.h'
      integer    num_children
      integer    i, sibling_in
      integer    temp_node(0:15)
      integer    temp_node_id
      integer    temp_parent(0:15)
      integer    local_seed
      integer    divisor
      integer    quotient
      integer    min_children
      integer    error
      integer    retval
C  MAX is defined in par_dfs.h   
      divisor = MAX((max_depth + 1)/max_children, 1)
C
      local_seed = parent(Seed)
      retval = My_rand(local_seed)
      num_children = MOD(retval  , max_children)
      quotient = parent(Depth)/divisor
      min_children = max_children - quotient - 1
      if (num_children .LT. min_children) then
          num_children = min_children
      endif
C
      call Copy_to_scratch(parent, temp_parent )
      local_seed = My_rand(local_seed + 
     +            My_rand(parent(Depth)))
      do 100 i = num_children- 1, 0, -1
      if (TREE_DEBUG) then
          print *, 'PTS Expand, i = ', i
      endif
          local_seed = My_rand(local_seed)
          call Initialize_node(temp_node, temp_parent,
     +                         i, local_seed)
          call Push(temp_node, local_stack)
          stack_size = stack_size + 1
 100   continue
C
       return
       end  
C
C *******************************************************************  
      subroutine Scatter( node_id, comm)
      include  'mpif.h'
      include  'node_stackf.h'
      include  'mainf.h'
C
      integer   node_id(0:15)
      integer   comm        
C
      integer  node(0:15)
      integer  q
      integer    status(MPI_STATUS_SIZE)
      integer    error
      integer    MAX_NODE_SIZE
C
      call MPI_COMM_SIZE(comm,  p, ierr)
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      if (my_rank .EQ. 0)  then
C  Don't use MPI_Scatter because we will need an extra   
C      buffer to avoid overwriting local_stack           
          do 100 q = p-1 , 1 , -1 
              call Pop(local_stack, node)
              call MPI_SEND(node, NodeSz,
     +            MPI_INTEGER, q,
     +             0, comm, ierr)
 100      continue
      else
C         MAX_NODE_SIZE = NODE_MEMBERS + max_depth - 1
          call MPI_RECV(node_id, NodeSz,
     +         MPI_INTEGER, 0, 0,
     +         comm,  status, ierr)
      endif
C
      return
      end
C
C
C *******************************************************************  
      subroutine Initialize(node, p, my_rank)
      integer   node(0:15)
      include  'mainf.h'
      include 'node_stackf.h'
C
      if (my_rank .NE. 0)  then
          call Push(node, local_stack)
      endif
C
      return
      end
