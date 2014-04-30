C  mainf.h - header file for main.f, and general use for global data
C
C  Global variables
       integer       retval
       integer       ierr
       integer       p
       integer       my_rank
       integer       io_comm
C
       integer       max_children
       integer       max_depth
       integer       cutoff_depth
       integer       max_work
        integer      StkSize
       common /Inpt/  max_children, 
     +               max_depth, 
     +               cutoff_depth, max_work,
     +               StkSize
C
       logical       STATS, TREE_DEBUG
       integer       StkMax
       integer       INFINITY
       integer       TRUE, FALSE
       parameter     (StkMax = 5000,
     +                INFINITY = 1000000,
     +                TRUE = 1, FALSE = 0)
       parameter (STATS = .TRUE., TREE_DEBUG = .FALSE.)
C
C functions
       integer Allocate_lists
       integer Allocate_type_arrays
       integer Calc_size
C
