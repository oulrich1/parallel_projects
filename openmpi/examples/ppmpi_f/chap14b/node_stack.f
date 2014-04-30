C  node_stack.f -- functions for manipulating tree nodes and the stack.  For
C      use in parallel tree search.
C 
C  See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
C      search.
C 
C *******************************************************************  
      subroutine Set_stack_scalars(count, stack)
      include 'node_stackf.h'
      include 'mainf.h'
      integer      count   
      integer      stack(0:StkSize)
      integer      subtract, i, Topptr
C
      stack(In_use) = count
      if (count .GE. NodeSz) then
             subtract = NodeSz
      else
             subtract = count
      endif
      stack(Top) = stack_list + count - subtract
      Topptr = stack(Top)
C
      if (TREE_DEBUG) then
         print *, 'Set Scalars, new top = : ', stack(Top),
     +            ' count = ', count, '1stval = ',
     +             stack(stack_list)
      endif
C zero out end of stack to remove any data prior to send
      do 100 i = stack(Top) + NodeSz, StkSize
           stack(i) = -1
 100  continue
      return
      end
C
C
C *******************************************************************  
      subroutine Copy_node(node1, node2, stack)
      include 'node_stackf.h'
      include 'mainf.h'
      integer node1 
      integer node2
      integer stack(0:StkSize)
      integer   i
C
      if (TREE_DEBUG) then
        print 10, 'Copy node ',node1, ' to node ', node2
        print 20, 'Num of elemets to copy = ',
     +              NodeSz
 10     format(A,I3,A,I3)
 20     format(A,I3)
      endif
      do 100 i = 0, NodeSz-1
           stack(node2 +i) = stack(node1 +i)
 100  continue
C
      return
      end
C
C
C *******************************************************************  
      subroutine Copy_to_scratch(node, scratch)
      integer  node(0:15) 
      integer  scratch(0:15)
C
      include 'node_stackf.h'
      integer   i
C
      do 100 i = 0, NodeSz - 1
           scratch(i) = node(i)
 100  continue
      return
      end
C
C
C *******************************************************************
C  Return 0 if successful, negative otherwise   
      integer function Allocate_lists( )
      include 'mainf.h'    
      include 'node_stackf.h'
C
      integer max_size_calc
      integer i
C
      max_size_calc = (max_depth+1)*(NODE_MEMBERS-1)
      max_size_calc = max_size_calc + max_depth *
     +                   (max_depth +1)/2
      max_size_calc = (max_children -1)*max_size_calc
      max_size_calc = max_size_calc + (NODE_MEMBERS - 1)
     +                    + max_depth
      max_size_calc = max_size_calc * NodeSz
      if (max_size_calc .GT. StkSize) then
         print *, '***ERROR: Increase StkSize in main.f',
     +            ' and node_stackf.h to: ',
     +             max_size_calc
         Allocate_lists = -1
      else
C initialize local stack values to zeroes initially
          do 100 i = 0, max_size_calc-1
              local_stack(i) = 0
 100      continue
          local_stack(Max_size) = max_size_calc
          local_stack(Top) = -1
          local_stack(In_use) = 0
C
          Allocate_lists = 0
      endif
C
      return
      end
C *******************************************************************  
      integer function Allocate_root(root_ptr)
      integer root_ptr
      include 'node_stackf.h'
      integer i
C
      local_stack(Top) = Stack_list
      root_ptr = local_stack(Top)
      Allocate_root = 0
C
      return
      end
C
C *******************************************************************  
      subroutine Initialize_node(node, parent,
     +           sibling_rank_in, seed_in)
      integer node(0:15)
      integer parent(0:15)
      integer     sibling_rank_in   
      integer     seed_in           
      integer  i, j
      include 'node_stackf.h'
C
      node(Depth) = parent(Depth) + 1
      node(Sibling_rank) = sibling_rank_in
      node(Seed) = seed_in
C 
      j = 0
      do 100 i = 0 , parent(Depth)-1
          node(Ancestor +i) = parent(Ancestor + i)
          j = j+1
 100  continue
      node(Ancestor + j) = parent(Sibling_rank)
       do 200 i = Ancestor + j + 1, NodeSz - 1
           node(i) = 0
 200   continue
C
      return
      end
C *******************************************************************  
      integer function Empty(stack)
      include 'node_stackf.h'
      include 'mainf.h'
      integer stack(0:StkSize)
C
      Topptr = stack(Top)
      if (stack(Top) .EQ. -1) then
          Empty =  TRUE
      else
          Empty  = FALSE
      endif
      if (TREE_DEBUG .AND. Empty .EQ. TRUE) then
          print *, 'In Empty, Stack Empty!'
      endif
      return
      end
C
C
C *******************************************************************  
C  Assumes Empty has already been called   
      subroutine Pop(stack, node)
      include 'mainf.h'
      include 'node_stackf.h'
      integer stack(0:StkSize)
      integer node(0:15)
      integer i
      integer Prev
C
      Topptr = stack(Top)
      do 100 i = 0, NodeSz -1
          node(i) = stack(Topptr + i)
          stack(Topptr + i) = 0
 100  continue
      if (stack(Top) .EQ. Stack_list) then
C we have popped last node, stack is empty
           stack(Top) = -1
      else
           Prev = stack(Top) - NodeSz
           stack(Top) = Prev
      endif
      stack(In_use) = stack(In_use) - NodeSz
C
      return
      end  
C
C
C *******************************************************************  
      subroutine Push(node,stack)
      include 'mainf.h'
      include 'node_stackf.h'
      integer  node(0:15)
      integer  stack(0:StkSize)
      integer i
C
      if (stack(Top) .EQ. -1) then
          stack(Top) = stack_list
      else
C move Top to end of last node on stack
          stack(Top) = stack(Top) +  NodeSz
      endif
      Topptr = stack(Top)
      do 100 i = 0, NodeSz-1
           stack(Topptr + i) = node(i)
 100  continue
      stack(In_use) = stack(In_use) + NodeSz
      return
      end
C
C
C *******************************************************************  
C  Assumes all processes have access to stdout   
      subroutine Print_node(nodeptr, stack, comm)
      include 'node_stackf.h'
      include 'mainf.h'
      include 'mpif.h'
      integer nodeptr
      integer stack(0:StkSize)
      integer comm, i
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      print *,'----------------------------------------------------'
      print 10,'Process ',my_rank, ' > Depth = ',
     +   stack(nodeptr + Depth),
     +   ' Sibling rank = ',
     +   stack(nodeptr +Sibling_rank),
     +  ', Seed = ',
     +    stack(nodeptr + Seed),', Size = ',
     +    NodeSz
 10   format(A, I3, A, I3, A, I3, A, I3, A, I3)
      print 20, 'Process ',my_rank,' > Ancestors = '
 20   format (A, I3, A)
      print *, (stack(nodeptr +Ancestor + i),
     +          i = 0, stack(nodeptr +Depth -1) )
      return
      end
C
C
C *******************************************************************  
      subroutine Print_stack(title, stack, comm)
      include 'node_stackf.h'
      include 'mainf.h'
      include 'mpif.h'
      character *20     title
      integer   stack(0:StkSize)
      integer   comm
C
      integer Prev
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      print *, '*****************************************************'
      print 10,'Process ',my_rank,' > ',title,
     +         '  Stack = '
 10   format (A, I3, A, A)
      Prev = stack(Top)
      do while (stack(Prev) .GE. Stack_list) 
          call Print_node(Prev, stack, comm)
C get Predecesor node
          Prev = Prev - stack(Prev + NodeSz)
      end do
      print *,'******************************************************'
      return
      end
C
C
C *******************************************************************  
      subroutine Print_stack_list(count, stack, comm)
      include 'node_stackf.h'
      include 'mainf.h'
      include 'mpif.h'
      integer     count   
      integer     stack(0:StkSize)
      integer     comm
C
      integer  i
      integer  s_ptr
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
      print 50, 'Process ',my_rank,' > count = ',count,
     +         ' Stack_list = '
 50   format(A, I3, A, I3, A)
      s_ptr = stack(stack_list)
      do 100 i = 0,  count-1 
          print 200, s_ptr
 200      format(I5)
          s_ptr = stack(stack_list + i)
 100  continue
C
      return
      end
C*******************************************************************
      integer function My_rand(seed)
      integer seed
C
      integer MODULUS, MULTIPLIER, INCREMENT
      parameter (MODULUS = 65537, MULTIPLIER = 16383,
     +              INCREMENT = 1001)

      My_rand = MOD( (MULTIPLIER*(seed) + INCREMENT) ,
     +                MODULUS)
      return
      end
