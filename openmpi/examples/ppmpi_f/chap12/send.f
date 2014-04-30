C  send.f -- wrapper for MPI_Send that takes timing information.
C      Should be compiled and linked with ./bcast.c
C 
C  Input: none
C  Output: Print_send_time prints total elapsed time spent in MPI_Send 
C 
C  See Chap 12, pp. 271 & ff in PPMPI.
C
      subroutine MPI_SEND(buffer, count, datatype,dest,tag,comm, ierr)
      include 'mpif.h'
      double precision send_time
      common send_time
      data send_time /0.0/
      integer buffer(*)
      integer count
      integer datatype
      integer dest
      integer tag
      integer comm
      integer ierr
C
      integer ierr1
      real start_time
      real finish_time
C
      start_time = MPI_WTIME(ierr1)
      call PMPI_SEND(buffer, count, datatype,
     +                  dest, tag, comm, ierr)
      finish_time = MPI_WTIME(ierr1)
      send_time = send_time + finish_time - start_time
C
      return
      end
C
      subroutine Print_send_time()
      include 'mpif.h'
      integer my_rank
      integer ierr
      double precision send_time
      common  send_time
C
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      print 100, 'Process ',my_rank, 
     +    ' > Total time in call MPI_SEND = ', 
     +    send_time, 'seconds'
 100  format (A, I3,  A, D14.8)
      return
      end
