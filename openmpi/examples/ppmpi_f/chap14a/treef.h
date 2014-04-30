C  treef.h -- header file for use in recursive and iterative dfs
C      programs
C 
C  MAX_NODES should be chosen so that it's at least as large as
C       1 - c^(d+1)
C      --------    
C       1 - c
C  where c = MAX_CHILDREN and d = max_depth
C
      integer MAX_CHILDREN 
      integer MAX_NODES 
      integer MAX_VALUE  
      integer FALSE 
      integer TRUE 
      parameter (MAX_CHILDREN = 3, MAX_NODES=1093, 
     +           MAX_VALUE = 1024, FALSE = 0, TRUE = 1)
C
C Following is a pseudo-structure for Tree_Node_t
C     integer  TREE_NODE_T(0:6)
      integer  value
      integer  parent
      integer  depth
      integer  child_count
      integer  child(0:2)
      parameter (value = 0, parent = 1, depth=2, 
     +           child_count = 3)
      data child /4,5,6/
C     
C     an array of tree_nodes_t's
      integer  tree(0:1092, 0:6)
      common /TreeBlk/  tree
      integer  tsize
      parameter (tsize = 1093 * 7)
C 
      integer MAX_STACK
      parameter (MAX_STACK = 20)
C  
C   a pseudo-structure for NODE_T
C     integer NODE_T(0:19)    
      integer size
      parameter (size = 0)
      integer list(20)
C
