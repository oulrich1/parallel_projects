C  node_stackf.h
C 
C  Basic definitions and declarations for stack and tree node manipulation
C
C  define struct for NODE_MEMBER - elements below
       integer NODE_MEMBERS
C  all nodes are the same size
       integer NodeSz
       parameter (NodeSz = 16)
       integer Depth, Sibling_rank, Seed 
       integer  Ancestors, Ancestor
       parameter (Depth = 0, Sibling_rank = 1, Seed = 2,
     +            Ancestors = 3, Ancestor = 3,
     +            NODE_MEMBERS = 4)
C
C typedef struct for Stack -- identifies elements
C stack TOP points to BEGINNING of LAST node on stack
        integer     max_size      
        integer     in_use
        integer     Top
        integer     stack_list
        parameter (max_size = 0, in_use = 1,
     +             Top =2, stack_list = 3)
        integer local_stack 
        integer scratch_node
        common /Stkcomm/ local_stack(0:5000),
     +                   scratch_node(0:15)
        integer Topptr
C 
C functions
        integer  Allocate_root
        integer  Empty
        integer  My_rand
        integer  Stack_int
C
C
