C  service_requests.f -- functions for handling requests for work
C     from other processes.
C 
C  See Chap 14, p. 332, in PPMPI. 
C
C *******************************************************************
C  Simply returns allocated size in Fortran
       integer function Allocate_type_arrays( )
       include 'mainf.h'
       include 'service_requestsf.h'
C
       allocated_size = cutoff_depth*(max_children - 1)/2
       Allocate_type_arrays = 0
       return
       end
C
C
C *******************************************************************  
      subroutine Service_requests( comm)
      integer comm          
C
      include 'mainf.h'
      include 'node_stackf.h'
      include 'service_requestsf.h'
      include 'queuef.h'
      integer destination
C
      do while (Work_requests_pending(comm) .NE. FALSE)
          destination = Get_dest(comm)
          if (Nodes_available() .NE. FALSE)  then
              call Split()
              call Send_work(destination, comm)
          else
              call Send_reject(destination, comm)
          endif
      end do 
C
      return
      end 
C
C
C *******************************************************************  
C  Count until at least two nodes above cutoff depth are found   
C      Should probably be combined with Split                    
      integer function Nodes_available()
C
      include 'node_stackf.h'
      include 'mainf.h'
      integer  node
      integer  Prev
      integer  count
C
      count = 0
      Prev = local_stack(Top)
C bottom of stack when Ptr = Stack_list
      do while ( (local_stack(Prev) .GE. Stack_list) .AND.
     +                         (count .LT. 2)) 
          if (local_stack(Prev + Depth) .LE. cutoff_depth) then
              count = count + 1
          end if
          Prev =  Prev - NodeSz
      end do
C
      if (count .GE. 2)then
          Nodes_available = TRUE
      else
          Nodes_available = FALSE
      endif
      return
      end   
C
C
C *******************************************************************  
C  Split builds a derived datatype that picks out the nodes to be
C      sent   
      subroutine Split()
      include 'mpif.h'
      include 'mainf.h'
      include 'node_stackf.h'
      include 'service_requestsf.h'
C
      integer     index
      integer     odd
      integer     node
      integer     i
      integer     Next
C
      node = Stack_list
      index = 0
      odd = 0 
      node_count = 0
      Topptr = local_stack(Top)
C 
      do while ( (local_stack(node) .NE. -1)
     +            .AND.   
     +           (node .LE. Topptr) )  
          if (local_stack(Depth + node) .LE. cutoff_depth) then
              if (odd .EQ. TRUE)  then
                  block_lengths(node_count) = NodeSz
                  displacements(node_count) = index
                  node_list(node_count) = node
                  node_count = node_count + 1
                  odd = 0     
              else
                  odd = 1     
              endif
          endif
          index = index +  NodeSz
          Next = node + NodeSz
          node = Next
      end do
C
      call MPI_TYPE_INDEXED(node_count, block_lengths,
     +     displacements, 
     +     MPI_INTEGER,
     +      send_stack_mpi_t, ierr)
      call MPI_TYPE_COMMIT( send_stack_mpi_t, ierr)
C
      return
      end   
C
C
C *******************************************************************  
      subroutine Send_work(destination, comm)
      integer       destination  
      integer       comm
C
      include 'mpif.h'
      include 'mainf.h'
      include 'queuef.h'
      include 'node_stackf.h'
      include 'statsf.h'
      include 'service_requestsf.h'
C
      if (TREE_DEBUG) then
          print *, '*** Send Work to ',destination,
     +             local_stack(stack_list+16),
     +             'Count = ',node_count*16
      endif
C This new type sends every other node in the list
      call MPI_SEND(local_stack(stack_list), 1,
     +                   send_stack_mpi_t, destination,
     +                   WORK_TAG, comm, ierr)
      call MPI_TYPE_FREE(send_stack_mpi_t, ierr)
      call Send_half_energy(destination, comm)
      if ( STATS) then
         call Incr_stat(work_sent)
      endif
C Now compress our local stack - delete nodes sent
      call Compress()
      return
      end   
C
C
C *******************************************************************
      subroutine  Compress()
      include 'mainf.h'
      include 'node_stackf.h'
      include 'service_requestsf.h'
      integer  compress_point, prev_compress                  
      integer  delete_node
      integer  save_node
      integer    new_in_use
      integer    i
      integer   Next
      integer   Topptr
      logical   update
C
      compress_point = node_list(0)
      prev_compress = compress_point
      new_in_use = displacements(0)
      update = .FALSE.
      if (TREE_DEBUG) then
          print 10,'Before Compress: top =',
     +             local_stack(Top),
     +            ' in use = ',local_stack(In_use),
     +            ' 1st comp pt =',
     +            compress_point
 10       format(A,I3,A,I3,A,I3)
      endif
C
      do 100 i = 0 , node_count -1
          delete_node = node_list(i)
          Next = delete_node + NodeSz
          if (Next  .LE. local_stack(Top))  then
      if (TREE_DEBUG) then
          print *, 'Compress:  comp =',compress_point,
     +             ' del =', delete_node
      endif
              save_node = Next
              call Copy_node(save_node, compress_point,
     +                            local_stack)
              new_in_use = new_in_use + NodeSz
              prev_compress = compress_point
              Next = compress_point + NodeSz
              compress_point = Next
              update = .TRUE.
          endif
 100  continue    
C
      if (update) then
          local_stack(Top) = prev_compress
          local_stack(In_use) = new_in_use
C  remove old data by filling in with zeros
          Topptr = local_stack(Top)
          do 200 i = Topptr + NodeSz,
     +               local_stack(StkSize)
                 local_stack(i) = -1
 200      continue
      endif
C
      if (TREE_DEBUG) then
          print *,'After Compress: top =',
     +            local_stack(Top),
     +            'in use = ',local_stack(In_use)
C          call Print_stack_list(local_stack(In_use),
C     +         local_stack, io_comm)
      endif
C
      return
      end   
C
C
C *******************************************************************  
      subroutine    Send_reject(destination, comm)
      integer       destination   
      integer       comm 
      integer       x
      include       'mpif.h'
      include       'mainf.h'
      include       'queuef.h'
      include       'statsf.h'       
C
      x = -1
C
      if (TREE_DEBUG) then
             print *, 'Send reject to ',destination
      endif
      call MPI_SEND( x, 1, MPI_INTEGER, destination,
     +            WORK_TAG, comm, ierr)
      if (STATS) then
         call Incr_stat(rejects_sent)
      endif
      return
      end
C
C
C *******************************************************************  
      subroutine Send_all_rejects(comm)
      integer comm 
      integer destination
      include 'mainf.h'
      include 'service_requestsf.h'
      include 'queuef.h'
      integer retval
C
      do while (Work_requests_pending(comm) .NE. FALSE)
          destination = Get_dest(comm)
          call Send_reject(destination, comm)
      end do
      return
      end
