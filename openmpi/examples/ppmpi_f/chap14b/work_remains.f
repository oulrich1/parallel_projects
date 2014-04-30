C  work_remains.f -- controls exit from main loop in parallel tree
C      search program.  Checks whether local stack is empty.  If not,
C      returns TRUE.  If it is, executes the following algorithm.
C 
C          Send "energy" to process 0 (see terminate.c, below).
C          while(TRUE) {
C              Send reject messages to all processes requesting work.
C              if (the search is complete) {
C                  Cancel outstanding work requests.
C                  return FALSE.
C              } else if (there's no outstanding request) {
C                  Generate a process rank to which to send a request.
C                  Send a request for work.
C              } else if (a reply has been received) {
C                  if (work has been received) return TRUE.
C              }
C          }
C 
C  See Chap 14, pp. 332 & ff, in PPMPI.
C
C *******************************************************************  
      integer function Work_remains( comm, p, my_rank)
      integer    comm          
C
      include      'work_remainsf.h'
      include      'mainf.h'
      include      'node_stackf.h'
      include      'statsf.h'
      integer      work_available
      integer      work_request_process
      integer      request_sent
      integer      cont
      integer      debug_count
C
      if (STATS) then
          call Start_time_rtn()
      endif
C
      cont = TRUE
      if (Empty(local_stack) .EQ. FALSE)  then
          if ( STATS) then
             call Finish_time(work_rem_time)
          endif
C
          retval = TRUE
      else  
          call Return_energy(comm)
          request_sent = FALSE
          do while (cont .EQ.TRUE  )
              call Send_all_rejects(comm)
              if (Search_complete(comm, p, my_rank) .EQ.
     +                TRUE)  then
C
                  if (STATS) then
                       call Finish_time(work_rem_time)
                  endif
                  if (request_sent .EQ. TRUE) then
                      call Cancel_request()
                  endif
                  retval =  FALSE
                  cont = FALSE
              else
                  if (request_sent .NE. TRUE)  then
                       work_request_process
     +                          = New_request(comm,p, my_rank)
                       call Send_request
     +                         (work_request_process, comm)
                       request_sent = TRUE
                  else
                       if (Reply_received(work_request_process,
     +                        work_available, comm)
     +                        .EQ. TRUE) then
                         if (work_available .EQ. TRUE)  then
                            if ( STATS ) then
                               call Finish_time(work_rem _time)
                            endif
                            retval = TRUE
                            cont = FALSE
                         else 
                            request_sent = FALSE
                         endif
                       endif
                   endif
C
              endif
C
          end do
      endif
      Work_remains = retval
      return
      end
C
C
C *******************************************************************  
      integer function New_request(comm, p, my_rank)
      integer comm
C
      include 'mainf.h'
      include 'work_remainsf.h'
      integer  rank
      integer  iseed
C
C
      iseed = my_rank+5
      rank = MOD(Random(iseed), p)
      do while (rank .EQ. my_rank)
         rank = MOD(Random(iseed) ,p)
      end do
C
      New_request = rank
      return 
      end   
C
C
C *******************************************************************  
      subroutine Send_request(work_request_process, comm)
      integer       work_request_process   
      integer       comm                   
C
      include 'work_remainsf.h'
      include 'mainf.h'
      include 'queuef.h'
      include 'node_stackf.h'
      include 'statsf.h'
      include 'mpif.h'
      integer     x
      data x /0/
C
      call MPI_SEND( x, 1, MPI_INTEGER, work_request_process,
     +               REQUEST_TAG, comm, ierr)
C
C  Post nonblocking receive
      call MPI_IRECV(local_stack(stack_list),
     +     StkSize,
     +     MPI_INTEGER, work_request_process, WORK_TAG, comm,
     +     posted_recv, ierr)
      if (TREE_DEBUG) then
         print *,'IRECV posted = ', 
     +           ' from = ',work_request_process,
     +           ' size = ',StkSize
         print *,'  '
      endif
      if (STATS) then
           call Incr_stat(requests_sent)
      endif
      return
      end
C
C
C *******************************************************************
      integer function Reply_received(work_request_process,
     +                work_available,  comm)
      integer      work_request_process   
      integer      work_available        
      integer      comm
C
      include 'work_remainsf.h'
      include 'mainf.h'
      include 'node_stackf.h'
      include 'statsf.h'
      include 'mpif.h'            
      logical      reply_rcvd
      integer      status(MPI_STATUS_SIZE)
      integer      count
C
      call MPI_TEST( posted_recv, reply_rcvd,  status, ierr)
C
      if (reply_rcvd)  then
      if (TREE_DEBUG) then
         print *, 'Reply Received. val =',
     +            local_stack(stack_list)
      endif
          if ( local_stack(stack_list) .EQ. -1)  then
              if (STATS) then
                 call Incr_stat(rejects_recd)
              endif
              work_available = FALSE
          else
              call MPI_GET_COUNT(status, MPI_INTEGER, count,ierr)
              call Set_stack_scalars(count, local_stack)
              call Recv_half_energy(work_request_process, comm)
              if (STATS) then
                  call Incr_stat(work_recd)
              endif
              work_available = TRUE
          endif
          Reply_received = TRUE
      else
          work_available = FALSE
          Reply_received = FALSE
      endif
      return
      end
C
C
C *******************************************************************  
C  Tree search has completed, but there is still an outstanding work   
C     request.  Try to cancel it.                                      
      subroutine Cancel_request()
      include 'mainf.h'
      include 'mpif.h'
      include 'work_remainsf.h'
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
C
C   Cancel not implemented in mpich!
      call MPI_CANCEL( posted_recv, ierr)
      call MPI_WAIT( posted_recv,  status, ierr)
      if (TREE_DEBUG) then
         print *,' Cancel = ',posted_recv
      endif
C
      return
      end
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
