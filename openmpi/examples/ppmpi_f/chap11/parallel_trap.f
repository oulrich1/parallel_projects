C  parallel_trap.f -- parallel trapezoidal rule with timing functions
C 
C  Input:
C      a: left endpoint
C      b: right endpoint
C      n: number of trapezoids
C  Output:
C      Integral
C      Elapsed time in seconds (excluding I/O).
C 
C  Note:  f(x) is hardwired.
C 
C  See Chap 11, pp. 254 & ff in PPMPI.
C 
C
      PROGRAM ParTrap
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
      integer p
      integer my_rank
      real    a
      real    b
      integer n
      real    h
      real    integral
      real    total 
      integer local_n
      real    local_a
      real    local_b
      integer io_comm
      integer i
      real    overhead
      real    start, finish
      integer retval
      integer ierr
      character *40 prompt
C function:
      real Trap 
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval = Cache_io_rank(MPI_COMM_WORLD, io_comm)
C
      total = 0.0
      prompt = 'Enter a, b, and n                       '
      retval = Cread(io_comm, 3, prompt, a,b, n)
C
C  Estimate overhead   
      overhead = 0.0
      do 100 i = 0 , 99  
          call MPI_BARRIER(MPI_COMM_WORLD, ierr )
          start = MPI_WTIME(ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr )
          finish = MPI_WTIME(ierr)
          overhead = overhead + (finish - start)
 100  continue
      overhead = overhead/100.0
C
      call MPI_BARRIER(MPI_COMM_WORLD, ierr )
      start = MPI_WTIME(ierr)
C
      h = (b-a)/n      
      local_n = n/p    
C
C  Length of each process' interval of
C  integration = local_n*h.  So my interval
C  starts at:   
      local_a = a + my_rank*local_n*h
      local_b = local_a + local_n*h
C
C  Call the serial trapezoidal function   
      integral = Trap(local_a, local_b, local_n, h)
C
C  Add up the integrals calculated by each process   
      call MPI_REDUCE( integral,  total, 1, MPI_REAL ,
     +     MPI_SUM, 0, MPI_COMM_WORLD, ierr )
C
      call MPI_BARRIER(MPI_COMM_WORLD, ierr )
      finish = MPI_WTIME(ierr)
      prompt = 'Our Estimate is                          '
      retval = Cprint(io_comm,1, 'Our estimate is        ',
     +                total)
      prompt = 'Elapsed time in seconds                 '
      retval = Cprint(io_comm,1, 'Elapsed time in seconds',
     +     (finish - start) - overhead)
C
      call MPI_FINALIZE(ierr)
      end
C
C
C ******************************************************************  
      real function Trap(local_a, local_b, local_n, h)
      real  local_a    
      real  local_b    
      integer    local_n    
      real  h          
C
      real integral     
      real x
      integer i
C function
      real f   
C
      integral = (f(local_a) + f(local_b))/2.0
      x = local_a
      do 100 i = 1, local_n-1
          x = x + h
          integral = integral + f(x)
 100  continue
      integral = integral*h
      Trap = integral
      return
      end
C
C
C ******************************************************************  
      real function f(x)
C  Calculate f(x).    
      real x
      f = x * x
      return
      end
C

