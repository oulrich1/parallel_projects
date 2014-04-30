C  queue.f -- functions for handling pending messages -- for use in
C      parallel tree search
C 
C  See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
C      search.
C
C *******************************************************************  
C  Only pending messages should be requests for work, rejects, or      
C      solutions.  Use stack as scratch area for receives.             
C      WARNING:  WILL OVERWRITE CONTENTS OF STACK!                     
      subroutine Clean_up_queues(comm)
      integer         comm
C
      include 'mpif.h'
      include 'mainf.h'
      include 'queuef.h'
      include 'node_stackf.h'
      integer         done
      logical         message_pending
      integer         stack(0:StkMax)
      integer         status(MPI_STATUS_SIZE)      
C
      done = FALSE
      do while (done .EQ. FALSE)
          call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG,
     +         comm,  message_pending, status, ierr)
          if (message_pending )  then
              call MPI_RECV(stack(stack_list),  StkSize,
     +            MPI_INTEGER, MPI_ANY_SOURCE, 
     +            MPI_ANY_TAG, comm,  status, ierr)
          else  
              done = TRUE
          endif
      end do
      return
      end
C
C
C *******************************************************************  
      integer function Work_requests_pending(comm)
      integer   comm
C
      include 'queuef.h'
      include 'mpif.h'
      include 'mainf.h'
      integer   status(MPI_STATUS_SIZE)
      logical   message_in_queue
C
      call MPI_IPROBE(MPI_ANY_SOURCE, REQUEST_TAG,
     +     comm,  message_in_queue, status, ierr)
C
      if (message_in_queue )  then
          Work_requests_pending = TRUE
      else  
          Work_requests_pending = FALSE
      endif
      return
      end
C
C
C *******************************************************************  
C  Only called when a message is pending   
      integer function Get_dest(comm)
      integer   comm
C
      include 'mpif.h'
      include 'queuef.h'   
      integer   node_request
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
C
      call MPI_RECV(node_request, 1, MPI_INTEGER,
     +     MPI_ANY_SOURCE, REQUEST_TAG,
     +      comm, status, ierr)
      Get_dest =  status(MPI_SOURCE)
      return
      end
C
