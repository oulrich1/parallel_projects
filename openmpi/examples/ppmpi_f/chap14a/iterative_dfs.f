C  iterative_dfs.f -- serial iterative depth-first search.  Generates
C      a random tree and uses iterative depth-first search to visit
C      the nodes in the tree and find the leaf node with minimum
C      value.
C 
C  Input:
C      max_depth:  the maximum depth of a node in the tree.
C 
C  Output:
C      1. The tree
C      2. A list of the nodes visited by depth-first search
C      3. The value of the minimum node.
C 
C  Algorithm:
C      1. Get max_depth
C      2. Set each entry in array storing tree to NULL
C      3. Recursive function Generate_tree uses rand function
C         to generate random tree.
C      4. Use iterative depth-first search to find leaf with
C         minimum value.
C      5. Print minimum value.
C 
C  Data structures:
C      The tree is a statically allocated 2-dimensional array
C  Each row of the array corresponds to a tree node.
C  Each tree node has a value -- a randomly
C  generated int with the property that the value assigned
C  to a parent is less than the value assigned to each of
C  its children.  Each tree node also stores the subscript of
C  its parent (-1 for root), its depth, the number of its
C  children, and subscripts of its children.  The type NODE_T is
C  the subscript of the node in the tree array.
C      The stack is just an array of NODE_T together with
C  an int indicating the number of elements in the array.
C 
C  Notes:
C      1.  MAX_NODES (the maximum number of nodes in a tree) should be 
C  chosen so that it's at least as large as
C       1 - c^(d+1)
C      --------    
C       1 - c
C  where c = MAX_CHILDREN and d = max_depth
C      2.  The array storing the tree and the current "best solution"
C  are global.
C 
C  See Chap. 14, pp. 324 & ff., in PPMPI.
C
C
C ******************************************************************  
      PROGRAM Iter
      INCLUDE 'treef.h'
      integer best_solution
C
      integer     max_depth
      integer     node_count
C
      print *, 'What''s the maximum depth of a node?' 
      read *,  max_depth
C      
      node_count = 0
      call Generate_tree(max_depth, node_count, -1, -1)
C
      call Print_tree( node_count)
      print *, "*** End of Tree ***"
C
      best_solution = 2*MAX_VALUE
      call Dfs_stack(0, best_solution)
      print *, 'The best solution is ', best_solution
C
      end
C
C ******************************************************************  
C Adds current value in *node_count to tree
      subroutine Generate_tree( max_depth, node_count,
     +                          parent_depth, parent_in)
      INCLUDE 'treef.h'
      integer     max_depth      
      integer     node_count     
      integer     parent_depth   
      integer     parent_in 
C
      integer i
      integer p
      data i /0/
      integer node 
      integer pow
C
      integer Random
      integer seed
      seed = 5
C
      node = node_count
      call Alloc_node(node)
      tree(node, depth) = parent_depth+1
      if (parent_in .EQ. -1) then
          tree(node,value) = MOD(Random(seed), MAX_VALUE)
      else         
         pow = 1
         do 100 p = 1, tree(node, depth)
              pow = pow * 2
 100     continue  
         tree(node,value) = MOD(Random(seed), (MAX_VALUE/ pow))
     +                    + tree(parent_in, value)
      endif
      tree(node, parent) = parent_in
C
      if (tree(node, depth) .NE. max_depth)  then
          do while ( (i .LT. MAX_CHILDREN) .AND.
     +            MOD(Random(seed),MAX_CHILDREN) .NE.0 )
              node_count = node_count + 1
              tree(node,child_count) = tree(node,child_count)+1
              tree(node,child(i)) = node_count
              call Generate_tree( max_depth, node_count, 
     +            tree(node, depth), node)
              i = i + 1
          end do  
      endif
      if (node .EQ. 0) then
         node_count = node_count + 1
      endif   
      return
      end
C
C
C ******************************************************************  
      subroutine Alloc_node(node)
      INCLUDE 'treef.h'
      integer     node   
      integer i
C
      do 100 i = 0,2 
          tree(node, child(i)) = -1
 100  continue
      tree (node, child_count) = 0
      return 
      end   
C
C
C ******************************************************************  
      subroutine Print_tree(node_count)
      INCLUDE 'treef.h'
      integer     node_count   
      integer node, i
C
      print 100,'node_count = ', node_count
 100  format (1x, A, I5)
      do 200 node = 0 , node_count-1  
          print *,'---------------------------------------------------'
          print 300,'Node = ',node, '  Value = ', tree(node,value)
 300      format(1x, A, I6, A, I8)
          print 400, 'Parent = ', tree(node,parent),' Depth = ', 
     +        tree(node,depth), ' Child_count = ',
     +        tree(node,child_count)
 400      format(1x, A, I5, A, I5, A, I5)
          print *, 'Children = '
          do 500 i = 0 , tree(node, child_count)-1
              print 600, ' ', tree(node, child(i))
 600          format (1x, A, i5)
 500      continue
 200  continue
      print *,'---------------------------------------------------'
      return
      end
C
C ******************************************************************  
C ******************************************************************  
C  Iterative depth-first search using a stack   
      subroutine Dfs_stack(root, best_solution) 
      INCLUDE   'treef.h'
      integer   root
      integer   node 
      integer   stack(0:20)
C
      integer best_solution
      integer val
C functions
      integer Solution
      integer Feasible
      integer Evaluate
      integer Empty
      integer Pop
C
C  Allocate and initialize stack   
      call INITialize( stack)
C
C  Expand root; push children onto stack   
      call Expand(root, stack)
C
      do while ( Empty(stack) .EQ. FALSE) 
          node = Pop( stack)
          if (Solution(node) .EQ. TRUE)  then
             val = Evaluate(node)
             if (val .LT. best_solution) then
                  best_solution = val
             endif
          else 
             if (Feasible(node,best_solution)
     +                 .EQ. TRUE)  then
                 call Expand(node, stack)
             endif
          endif
      end do
C
      return
      end  
C
C
C ******************************************************************  
      subroutine INITialize(stack) 
C
      integer stack(0:20)
C
C     element 0 determines stack size
      stack(0) = 0
C
      return
      end   
C
C
C ******************************************************************  
      subroutine Expand(node, stack)
      INCLUDE 'treef.h'
      integer node     
      integer stack(0:20)   
C
      integer i
      integer num_children
C
      num_children = tree(node, child_count)
      i = num_children - 1
      do while ( i .GE. 0  )
          call Push( tree(node,child(i)), stack)
          i = i - 1
      end do
C
      return
      end   
C
C
C ******************************************************************  
      integer function Empty(stack) 
C
      integer stack(0:20)
C
C     element 0 = stack size.   true = 1, false = 0
      if (stack(0) .EQ. 0) then
          Empty = 1
      else
          Empty = 0
      endif
C
      return
      end
C
C
C ******************************************************************  
      subroutine Push( node,  stack) 
C
C   a psuedo-structure for NODE_T
      integer NODE_T(0:19)    
      integer size
      integer list(20)
      parameter (size = 0)
C
      integer node 
      integer stack(0:20)
C
      stack(stack(size)+1) = node
      stack(size) = stack(size) + 1
C
      return
      end   
C
C
C ******************************************************************  
      integer function Pop( stack)
C
      integer stack(0:20)
      integer node 
      integer size
      integer list(20)
      parameter (size = 0)  
C
      stack(size) = stack(size) - 1
      Pop = stack( stack(size) + 1)
      return
      end   
C
C
C ******************************************************************  
C  A solution node has no children   
      integer function Solution(node)
      INCLUDE 'treef.h'     
      integer     node   
C  1 = True, 0 = False
C
      if ( tree(node, child_count) .EQ. 0) then
          Solution = 1
      else
          Solution = 0
      endif   
C
      return 
      end
C
C ******************************************************************  
      integer  function Evaluate(node)
      INCLUDE  'treef.h'
      integer     node    
C
      Evaluate = tree(node,value)
      print *, "Evaluate = ", Evaluate
      return
      end   
C
C
C ******************************************************************  
      integer function Feasible(node, best_solution) 
      INCLUDE 'treef.h'
      integer     node 
      integer     best_solution
C
C  True = 1, False = 0          
      if (tree(node,value) .LT. best_solution) then
          Feasible = 1
      else
          Feasible = 0
      endif
C
      return
      end
C
C**********************************************************
C
      integer function Random(seed)
C
      integer seed
      integer c1, c2
      parameter (c1 = 19423, c2 = 3011)
      integer oldseed
      save oldseed
      data oldseed /0/
C
      if (oldseed .EQ. 0) then
          oldseed = seed
      endif
C
      oldseed = mod(c1*oldseed, c2)
C
      Random = oldseed 
C
      return
      end
