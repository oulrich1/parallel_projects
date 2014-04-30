C  par_dfs.f -- local depth-first search function.  Performs iterative
C      depth-first search on tree until local_stack is exhausted or
C      it has expanded max_work nodes.
C 
C  See Chap 14, pp. 330 & ff., in PPMPI.
C 
C
C *******************************************************************  
      subroutine Par_dfs(comm, p, my_rank)
      integer comm
C
      include 'mainf.h'
      include 'par_dfsf.h'
      include 'node_stackf.h'
      include 'statsf.h'
      include 'solutionf.h'
C
      integer     temp_sol
      integer     count, topptr
      integer     node(0:15)
C
C  Search local subtree for a while   
      count = 0
      do while (Empty(local_stack) .NE. TRUE .AND.
     +         (count .LT. max_work))
      if (TREE_DEBUG) then
          topptr = local_stack(Top)
          print *, my_rank, 'PAR_DFS Pop node',
     +    '.Top = ',local_stack(Top) 
      endif
          call Pop(local_stack, node)
          if (STATS) then
              call Incr_stat(nodes_expanded)
          endif
          if (Solution(node) .NE. FALSE )  then
              temp_sol = Evaluate(node)
      if (TREE_DEBUG) then
              print *, my_rank, ' In PAR_DFS Soln True:
     + Temp sol(0) = ', temp_sol
      endif
              if (temp_sol .LT.
     +                Best_solution_funct(comm)) then
                  call Local_solution_update(temp_sol, node)
                  call Bcast_solution(comm)
              endif
          else
      if (TREE_DEBUG) then
              print *,my_rank,'In PAR_DFS Soln False. Top=',
     +              local_stack(Top)
      endif
              if (Feasible(node, comm) .NE. FALSE) then
                   call Expand(node, p, my_rank)
      if (TREE_DEBUG) then
                   print *,my_rank,
     +             'Return from PAR_DFS Expand'
      endif
              endif
          endif
          count = count + 1
      end do
C
      return
      end
C
C
C ******************************************************************  
      integer function Solution(node)
      integer  node(0:15)
C
      include 'mainf.h'
      include 'node_stackf.h'
      if (node(Depth) .EQ. max_depth) then
          Solution = TRUE
      else
          Solution = FALSE
      endif
      return
      end
C
C
C ******************************************************************  
      integer function Evaluate(node)
      integer node(0:15)   
C
      include 'node_stackf.h'
C
      Evaluate = node(Seed)
      return
      end
C
C
C ******************************************************************  
      integer function Feasible(node, comm)
      integer node(0:15)
      integer comm
C
      Feasible = 1
      return
      end
C
C ******************************************************************  
      subroutine Expand(parent, p, my_rank)
      integer  parent(0:15)   
C
      include 'node_stackf.h'
      include 'par_dfsf.h'
      include 'mainf.h'
      include 'mpif.h'
C
      integer    num_children
      integer    i
      integer    temp_node(0:15)
      integer    temp_node_id
      integer    temp_parent(0:15)
      integer    local_seed
      integer    divisor
      save       divisor
      integer    quotient
      integer    min_children
      integer    error
      data       divisor /0/
      save       divisor
C
      if (divisor .EQ. 0) then
          divisor = MAXAB((max_depth + 1)/max_children, 1)
      endif
C
      local_seed = parent(Seed)
      num_children = MOD( My_rand(local_seed) ,max_children)
      quotient = parent(Depth)/divisor
      min_children = max_children - quotient - 1
      if (num_children .LT. min_children) then
          num_children = min_children
      endif
C
      call Copy_to_scratch(parent, temp_parent)
      local_seed = My_rand(local_seed +
     +               My_rand(parent(Depth)))
      i = num_children - 1
      do while (i .GE. 0)
          local_seed = My_rand(local_seed)
      if (TREE_DEBUG) then
          print *, my_rank, 'DFS Expand - i = ',i, local_seed
      endif
          call Initialize_node(temp_node,temp_parent,i,local_seed)
          call Push(temp_node, local_stack)
          i = i - 1
      end do
      return
      end
C   
C*****************************************************************
      integer function MAXAB(A, B)
      integer A, B
C
      if (A .GE. B) then
         MAXAB = A
      else
         MAXAB = B
      endif
      return
      end
