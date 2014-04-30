C  comm_time.f  
C  Version 2:  Process 0 starts off the ring pass.
C  Time communication around a ring of processes.
C      Guaranteed to have bugs.
C  
C  Input: None (see notes).
C 
C  Output:  Average, minimum, and maximum time for messages 
C     of varying sizes to be forwarded around a ring of 
C     processes.
C 
C  Algorithm:
C     1.  Allocate and initialize storage for messages 
C         and communication times
C     2.  Compute ranks of neighbors in ring.
C     3.  Foreach message size
C     3b.     Foreach test
C     3a.         Start clock
C     3c.         Send message around loop
C     3d.         Add elapsed time to running sum
C     3e.         Update max/min elapsed time
C     4.  Print times.
C 
C  Functions:
C     Initialize:  Allocate and initialize arrays
C     Print_results:  Send results to I/O process
C         and print.
C 
C  Notes:  
C     1. Due to difficulties some MPI implementations 
C        have with input, the number of tests, the max 
C        message size, the min message size, and the size 
C        increment are hardwired.
C     2. We assume that the size increment evenly divides
C        the difference max_size - min_size
C     3. Link with cio_mod1.o.
C 
C  See Chap 9, pp. 192 & ff and pp. 202 & ff in PPMPI.

      PROGRAM Comm2
      INCLUDE 'ciof.h'
      INCLUDE 'mpif.h'
      integer         test_count
      integer         max_size
      integer         min_size   
      integer         size_incr  
      parameter (max_size = 1000, min_size = 0,
     +           size_incr = 500, test_count = 2)
C      msg. sizes     
      real      x(0:999)
      double precision     times(0:999)
      double precision     max_times(0:999)
      double precision     min_times(0:999)
      integer  time_array_order
C      arrays.        
      double precision      start
      double precision      elapsed
      integer   i, test, size       
      integer   p, my_rank, source, dest
      integer   io_comm
      integer   status(MPI_STATUS_SIZE)
      integer   ierr, retval
C
      IO_KEY = MPI_KEYVAL_INVALID    
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval = Cache_io_rank(MPI_COMM_WORLD, io_comm)
C
      print 10, 'Proc: ', my_rank, ' P = ',p ,
     +                ' Before call Initialize   '
      call Initialize(max_size, min_size, size_incr, my_rank,
     +      x,  times,  max_times,  min_times,
     +      time_array_order)
C
      source = MOD((my_rank -  1) , p)
      dest = MOD((my_rank + 1) , p)
C
C  For each message size, find average circuit time   
C      Loop var size = message size                   
C      Loop var i = index into arrays for timings 
      size = min_size 
      i = 0    
      do while (size .LE. max_size)
          print 10, 'Proc: ', my_rank, ' P = ',p ,
     +                ' Before IF MYRANK   '
          if (my_rank .EQ. 0)  then
             times(i) =0.0
             max_times(i) = 0.0
             min_times(i) = 1000000.0
             do 100 test = 0 ,test_count-1  
                 start = MPI_WTIME(ierr)
                  call MPI_SEND(x, size, MPI_REAL , dest,
     +              0, MPI_COMM_WORLD, ierr )
                  call MPI_RECV(x, size, MPI_REAL , source,
     +              0, MPI_COMM_WORLD,  status, ierr )

                 elapsed = MPI_WTIME(ierr) - start
                 times(i) = times(i) + elapsed
                 if (elapsed .GT. max_times(i)) then
                     max_times(i) = elapsed
                 endif
                 if (elapsed .LT. min_times(i)) then
                     min_times(i) = elapsed
                 endif
 100          continue 
          else    
              do 200 test = 0 ,test_count-1      
                  call MPI_RECV(x, size, MPI_REAL , source, 
     +                  0, MPI_COMM_WORLD,  status, ierr )
                  call MPI_SEND(x, size, MPI_REAL , dest, 
     +                 0, MPI_COMM_WORLD, ierr )
 200          continue             
         endif
         size = size + size_incr
         i = i + 1
      end do    
C
      call Print_results(io_comm, my_rank, min_size, max_size,
     +     size_incr, time_array_order, test_count, times,
     +     max_times, min_times)
C
 10   format(A, I3, A, I3, A)
      call MPI_FINALIZE(ierr)
      end
C
C ******************************************************  
      subroutine Initialize(max_size, min_size,  size_incr,
     + my_rank, x, times,
     + max_times, min_times,order) 
      integer max_size
      integer min_size
      integer size_incr
      integer my_rank
      real    x(0:999)
      double precision    times(0:999)
      double precision    max_times(0:999)
      double precision    min_times(0:999)
      integer order
      integer i
C
      order = (max_size - min_size)/size_incr
C
C  Initialize buffer -- why this?   
      do 100 i = 0, max_size-1
          x(i) = my_rank
 100  continue  
      return
      end
C
C
C ******************************************************  
C  Send results from process 0 in MPI_COMM_WORLD to       
C  I/O process in io_comm, which prints the results.      
      subroutine Print_results(  io_comm,  my_rank,
     +      min_size,   max_size,   size_incr,
     +     time_array_order, test_count,   times,
     +     max_times,  min_times) 
      integer io_comm
      integer my_rank
      integer min_size, max_size
      integer size_incr
      integer time_array_order
      integer test_count
      double precision    times(0:999)
      double precision    max_times(0:999), min_times(0:999)
C
      include 'mpif.h'
      include 'ciof.h'
      integer i
      integer size
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      integer io_process
      integer io_rank
      integer retval
C
      retval= Get_io_rank(io_comm,  io_process)
      call MPI_COMM_RANK(io_comm,  io_rank, ierr)
C
      if (my_rank .EQ. 0)  then
          call MPI_SEND(times, time_array_order,
     +         MPI_DOUBLE_PRECISION,
     +         io_rank, 0, io_comm, ierr)
          call MPI_SEND(max_times, time_array_order,
     +         MPI_DOUBLE_PRECISION,
     +         io_process, 0, io_comm, ierr)
          call MPI_SEND(min_times, time_array_order,
     +         MPI_DOUBLE_PRECISION,
     +         io_process, 0, io_comm, ierr)
      endif
      if (io_rank .EQ. io_process)  then
          call MPI_RECV(times, time_array_order,
     +         MPI_DOUBLE_PRECISION,
     +         MPI_ANY_SOURCE, 0, io_comm,  status, ierr)
          call MPI_RECV(max_times, time_array_order,
     +         MPI_DOUBLE_PRECISION,
     +         MPI_ANY_SOURCE, 0, io_comm,  status, ierr)
          call MPI_RECV(min_times, time_array_order,
     +         MPI_DOUBLE_PRECISION,
     +         MPI_ANY_SOURCE, 0, io_comm,  status, ierr)
C
          print *,'Message size (floats):  '
          do 200 size =min_size, max_size-1, size_incr
              print *,  size
 200      continue
C
          print *,'Avg circuit time (ms):  '
          print *,((1000.0*times(i)/test_count),
     +             i = 0, time_array_order-1)
C
          print *,'Max circuit time (ms):  '
          print *,(1000.0*max_times(i), i = 0, 
     +                  time_array_order-1)
C
          print *,'Min circuit time (ms):  '
          print *, (1000.0*min_times(i),
     +                  i = 0, time_array_order-1)
      endif
C
      return
      end

