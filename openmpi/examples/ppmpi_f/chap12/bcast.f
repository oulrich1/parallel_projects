C  bcast.f -- linear loop broadcast, illustrates use of MPI's
C      profiling interface.
C 
C  Input: none
C  Output:
C      Total time spent by each process in MPI_Send
C      Result of broadcast -- should be "1" printed by each
C 
C  Note:  Should be linked with ./send.c and the MPI profiling
C      library.
C 
C  See Chap 12, pp. 271 & ff in PPMPI.
C
C
      PROGRAM BCast
      INCLUDE 'mpif.h'
      integer           p
      integer           my_rank
      integer           root
      integer           comm
      integer           proc
      integer           status(MPI_STATUS_SIZE)
      integer           ierr
      integer           size
      integer           tag 
      data   tag, size, root/ 0, 1, 0/
C
      integer          x
      integer          y
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  comm, ierr )
C
      if (my_rank .EQ. root) then
          y = 1
      else
          y = 0
      endif
      x =  y
C
      if (my_rank .EQ. root)  then
          do 100 proc = 0   ,p-1  
              call MPI_SEND(x, size, MPI_INTEGER, 
     +                      proc, tag, comm, ierr)
 100      continue
      endif
C
      call MPI_RECV(x, size, MPI_INTEGER, root, 
     +              tag, comm,  status, ierr)
C
      call Print_send_time()
      print 200 ,'Process ',my_rank, ' > x = ',  x
 200  format(A, I3, A, I3)    
      call MPI_FINALIZE(ierr)
      end
