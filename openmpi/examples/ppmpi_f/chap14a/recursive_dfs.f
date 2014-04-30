C  recursive_dfs.f -- serial recursive depth-first search.  Generates
C      a random tree and uses recursive depth-first search to visit
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
C      4. Use recursive depth-first search to find leaf with
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
C 
C  Notes:
C      1.  Fortran implementation must allow recursive subroutine
C  calls.
C      2.  MAX_NODES (the maximum number of nodes in a tree) should be 
C  chosen so that it's at least as large as
C       1 - c^(d+1)
C      --------    
C       1 - c
C  where c = MAX_CHILDREN and d = max_depth
C      3.  The array storing the tree, and the current "best solution"
C  are global.
C 
C  See Chap. 14, pp. 324 & ff., in PPMPI.
C
      PROGRAM Recrsv
      INCLUDE 'treef.h'
      data  tree / tsize*0/
C
      integer     best_solution
      integer     max_depth
      integer     node_count
C
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
      call Dfs_recursive(0, best_solution)
      print *, 'The best solution is ', best_solution
C
      end
C ******************************************************************  
C Adds current value in *node_count to tree
      subroutine Generate_tree( max_depth, node_count,
     +                          parent_depth, parent_in)
      INCLUDE 'treef.h'
      integer max_depth
      integer node_count
      integer parent_depth
      integer parent_in
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
C ******************************************************************
      subroutine Alloc_node(node)
      INCLUDE   'treef.h'
C
      integer     node   
      integer i
C
      do 100 i = 0,MAX_CHILDREN -1 
          tree(node, child(i)) = -1
 100  continue
      tree (node, child_count) = 0
      return 
      end   
C
C ******************************************************************  
      subroutine Print_tree(node_count)
      INCLUDE 'treef.h'
      integer  node_count
      integer  node
      integer  i
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
      subroutine Dfs_recursive(node, best_solution)
      INCLUDE   'treef.h'  
      integer    node
      integer    best_solution
      integer    num_children
      integer    child_list(0:MAX_CHILDREN-1)
      integer    i
      integer    val
C functions
      integer    Feasible
      integer    Solution
      integer    Evaluate
C
      print *,'Visiting Node ', node
C
      call Expand(node, child_list,  num_children)
C
      do 100 i = 0 , num_children-1
           print *,'Checking node ',node
           if (Solution(child_list(i)) .EQ. TRUE)  then
              val = Evaluate(child_list(i))
              if (val .LT. best_solution) then
                  best_solution =  val
              endif
          else
              if (Feasible(child_list(i),best_solution)
     +                  .EQ. TRUE)  then
                  call Dfs_recursive(child_list(i), best_solution)
              endif
          endif
 100  continue
      return
      end  
C
C ******************************************************************  
      subroutine Expand(node, child_list, num_children)
      INCLUDE 'treef.h'
      integer node     
      integer child_list(0:MAX_CHILDREN)
      integer num_children

      integer i 
C
      num_children = tree(node, child_count)
      do 100 i = 0, num_children-1
        child_list(i) = tree(node,child(i))
 100  continue 
C
      return
      end   
C
C ******************************************************************  
C  A solution node has no children   
      integer function Solution(node)
      INCLUDE   'treef.h'
      integer     node   
C
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
      integer     node    
      include     'treef.h'
C
      Evaluate = tree(node,value)
      print *, "Evaluate = ", Evaluate
      return
      end   
C
C ******************************************************************  
      integer function Feasible(node, best_solution) 
      integer     node 
      integer     best_solution
      INCLUDE   'treef.h'
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
