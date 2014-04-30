C  main.f
C 
C  Main program that calls Par_tree_search.
C 
C  Input:
C      1. max_depth: the maximum depth of nodes in the tree.
C      2. cutoff_depth:  the maximum depth of nodes that can be sent
C         to other processes in dynamic load balancing.
C      3. max_work:  the maximum number of nodes expanded in a
C         single call to Par_dfs.
C      4. max_children:  the maximum number of children that a
C         tree node can have.
C 
C  Output:
C      1. If STATS has been defined, statistics on the execution
C         of the program.  See stats.f
C 
C  Algorithm:
C      1. Start up MPI and get input.
C      2. Call various setup functions.
C      3. Process 0:  initialize root of tree
C      4. Call Par_tree_search
C      5. Print stats
C      6. Clean up message queues
C      7. Free storage
C 
C  Global Variables:
C      max_work: maximum number of nodes expanded in a single call
C          Par_dfs
C      cutoff_depth: maximum depth of nodes redistributed in load
C          balancing
C      max_children: maximum number of children of a node
C      max_depth: maximum depth of a node
C      io_comm: communicator for I/O -- duplicate of MPI_COMM_WORLD
C      p: number of processes in MPI_COMM_WORLD
C      my_rank: process rank in MPI_COMM_WORLD
C 
C  See Chap 14, pp. 328 & ff., in PPMPI.
C 
C
C *******************************************************************
      PROGRAM MainTree
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
      INCLUDE 'mainf.h'
      INCLUDE 'par_tree_searchf.h'
      INCLUDE 'node_stackf.h'
      INCLUDE 'solutionf.h'
      INCLUDE 'statsf.h'
      INCLUDE 'queuef.h'
      INCLUDE 'terminatef.h'
C
      integer       error
      integer       root
      integer       i
C
      IO_KEY = MPI_KEYVAL_INVALID
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval = Cache_io_rank(MPI_COMM_WORLD, io_comm )
C
C initialize our common data
      if (STATS) then
         call Set_stats_to_zero(statsi, statsd)
      endif
      returned_energy(0) = 0
      returned_energy(1) = 1
      ONE(0) = 1
      ONE(1) = 1
      StkSize = StkMax
C
C      call Cread(io_comm, 4,
C     +     'Enter max depth, cutoff depth, max work, & max children',
C     +      max_depth,  cutoff_depth,  max_work, max_children)
      max_depth = 10
      cutoff_depth = max_depth
      max_children = 4
      max_work = 3
      if (my_rank .EQ. 0) then
         print *,'max_depth, cutoff_depth, max_work, max_children: '
         print 10, max_depth, cutoff_depth, max_work, max_children
10       format(I3, I3, I3, I3)
      endif
      call Setup_term_detect(p, my_rank)
C
      error = Initialize_soln( )
      retval = Cerror_test(io_comm, 'Initialize_soln', error)
C
      error = Allocate_lists( )
      retval = Cerror_test(io_comm, 'Allocate_lists', error)
C
      error = Allocate_type_arrays()
      retval = Cerror_test(io_comm, 'Allocate_type_arrays', error)
C
C use -1 for "null"
      if (my_rank .EQ. 0)  then
          call Get_root( root)
      else  
          root = -1
      endif
C
      if (TREE_DEBUG) then
          print *, '***Begin Parallel Tree Search***'
      endif
      call Par_tree_search(root, io_comm)
C
      if (STATS) then
           call Print_stats(io_comm)
      endif
C 
      call Clean_up_queues(MPI_COMM_WORLD)
      call MPI_FINALIZE(ierr)
C
      end
C
C
C *******************************************************************  
      subroutine Get_root(root_ptr)
      integer  root_ptr
      integer  root(0:15)
      include 'node_stackf.h'
      include 'mainf.h'
C
      call Allocate_root(root_ptr)
      root( Depth) = 0
      root(Sibling_rank) = 0
      root(Seed) = 1
C
      call Push(root, local_stack)
C
      return
      end
