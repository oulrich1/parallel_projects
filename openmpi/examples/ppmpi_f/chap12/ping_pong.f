C  ping_pong.f -- two-process ping-pong -- send from 0 to 1 and send back
C      from 1 to 0
C 
C  Input: none
C  Output: time elapsed for each ping-pong
C 
C  Notes:
C      1.  Size of message is MAX_ORDER floats.
C      2.  Number of ping-pongs is MAX.
C 
C  See Chap 12, pp. 267 & ff. in PPMPI.
C 
C
      PROGRAM PingPong
      INCLUDE 'mpif.h'
      integer    p
      integer    my_rank
      integer    test
      integer    min_size  
      integer    max_size 
      integer    incr 
      parameter  (min_size = 0, max_size = 16, incr = 8)
      real       x(0:999)
      integer    size
      integer    pass
      integer    status(MPI_STATUS_SIZE)
      integer    ierr
      integer    i
      real       wtime_overhead
      real       start, finish
      real       raw_time
      integer    comm
      integer    MAX_ORDER, MAX
      parameter  (MAX_ORDER = 1000, MAX = 2)
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  comm, ierr )
C
      wtime_overhead = 0.0
      do 100 i = 0, 99
          start = MPI_WTIME(ierr)
          finish = MPI_WTIME(ierr)
          wtime_overhead = wtime_overhead + (start - finish)
 100  continue
      wtime_overhead = wtime_overhead/100.0
C
      if (my_rank .EQ. 0)  then
          test = 0
          size = min_size
          do while (size .LE. max_size)  
              do 200 pass = 0 , MAX -1
                  call MPI_BARRIER(comm, ierr)
                  start = MPI_WTIME(ierr)
                  call MPI_SEND(x, size, MPI_REAL , 1, 0, comm,ierr)
                  call MPI_RECV(x, size, MPI_REAL , 1, 0, comm,
     +                   status, ierr)
                  finish = MPI_WTIME()
                  raw_time = finish - start - wtime_overhead
                  print *, 'Size = ',  size,' Raw Time = ', raw_time
 200          continue
              size = size + incr
              test = test + 1
          end do
      else
          test = 0
          size = min_size
          do while (size .LE. max_size)  
              do 300 pass = 0 ,MAX - 1
  		  call MPI_BARRIER(comm, ierr)
                  call MPI_RECV(x, size, MPI_REAL , 0, 0, comm,
     +		     status, ierr)
                  call MPI_SEND(x, size, MPI_REAL , 0, 0, comm, ierr)
 300          continue
              size = size + incr
              test = test + 1
          end do
      endif
C
C
      call MPI_FINALIZE(ierr)
      end
