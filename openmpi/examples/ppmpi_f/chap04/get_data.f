C    get_data.f -- Parallel Trapezoidal Rule -- Uses basic 
C       Get_data function for input.
C
C    Input:
C       a, b: limits of integration.
C       n: number of trapezoids.
C    Output:  Estimate of the integral from a to b of f(x)
C       using the trapezoidal rule and n trapezoids.
C   
C    Algorithm:
C       1.  Each process calculates 'its' interval of
C           integration.
C       2.  Each process estimates the integral of f(x)
C           over its interval using the trapezoidal rule.
C       3a. Each process .NE. 0 sends its integral to 0.
C       3b. Process 0 sums the calculations received from
C           the individual processes and prints the result.
C   
C  Notes: 
C     1.  f(x) is hardwired.
C     2.  Assumes number of processes (p) evenly divides
C         number of trapezoids (n)
C   
C    See Chap. 4, pp. 60 & ff in PPMPI.
C
C
      PROGRAM getdat
      INCLUDE 'mpif.h'
      integer    my_rank     
      integer    p           
      real       a           
      real       b      
      integer    n           
      real       h           
      real       local_a     
      real       local_b     
      integer    local_n     
C  my calculation              
      real       integral    
      real       total       
      integer    source      
      integer    dest    
      integer    tag 
      integer    status(MPI_STATUS_SIZE)
      integer    ierr
      data dest,tag /0, 0/
C
C  function definitions
      real Trap      
      real f
C
C  Let the system do what it needs to start up MPI   
      call MPI_INIT( ierr)
C
C  Get my process rank   
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
C  Find out how many processes are being used   
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
C
      call Get_data( a,  b,  n, my_rank, p)
C
      h = (b-a)/n      
      local_n = n/p    
C
C  Length of each process' interval of 
C  integration = local_n*h.  So my interval starts at: 
      local_a = a + my_rank*local_n*h
      local_b = local_a + local_n*h
      integral = Trap(local_a, local_b, local_n, h)
C
C  Add up the integrals calculated by each process   
      if ( my_rank .EQ. 0)   then
          total = integral
          do 100 source = 1, p-1 
              call MPI_RECV( integral, 1, MPI_REAL , source, tag,
     +                       MPI_COMM_WORLD,  status, ierr )
	      total = total + integral
 100      continue
      else  
          call MPI_SEND( integral, 1, MPI_REAL , dest,
     +                   tag, MPI_COMM_WORLD, ierr )
      endif
C
C  Print the result   
      if  (my_rank .EQ. 0)   then
          print *, 'With n = ', n, ' trapezoids, our estimate'
          print *, 'of the integral from', a ,' to ',b,' = ',total
      endif
C
C  Shut down MPI   
      call MPI_FINALIZE(ierr)
      end
C
C
C ******************************************************************  
C  Function Get_data
C    Reads in the user input a, b, and n.
C    Input parameters:
C        1.  integer my_rank:  rank of current process.
C        2.  integer p:  number of processes.
C    Output parameters:
C        1.  float* a:  pointer to left endpoint a.
C        2.  float* b:  pointer to right endpoint b.
C        3.  int* n:  pointer to number of trapezoids.
C    Algorithm:
C        1.  Process 0 prompts user for input and
C            reads in the values.
C        2.  Process 0 sends input values to other
C            processes.
C   
      subroutine Get_data(a, b, n, my_rank, p)
      real  a     
      real  b 
      integer    n    
      integer    my_rank   
      integer    p
      INCLUDE 'mpif.h'
C
      integer source       
      integer dest            
      integer tag
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      data source /0/
C
      if  (my_rank .EQ. 0)  then
         print *, 'Enter a, b and n'
         read *, a, b, n
C
C
          do 100 dest = 1 ,  p-1  
	      tag = 0
              call MPI_SEND(a, 1, MPI_REAL , dest, tag,
     +                 MPI_COMM_WORLD, ierr )
              tag = 1
              call MPI_SEND(b, 1, MPI_REAL , dest, tag,
     +                 MPI_COMM_WORLD, ierr )
	      tag = 2
              call MPI_SEND(n, 1, MPI_INTEGER, dest,
     +                 tag, MPI_COMM_WORLD, ierr )
 100      continue
      else  
          tag = 0
          call MPI_RECV(a, 1, MPI_REAL , source, tag,
     +          MPI_COMM_WORLD,  status, ierr )
          tag = 1
          call MPI_RECV(b, 1, MPI_REAL , source, tag,
     +          MPI_COMM_WORLD,  status, ierr )
          tag = 2
          call MPI_RECV(n, 1, MPI_INTEGER, source, tag,
     +              MPI_COMM_WORLD,  status, ierr )
      endif
      return
      end
C
C
C ******************************************************************  
      real function Trap(local_a, local_b, local_n, h)
      real  local_a    
      real  local_b    
      integer local_n    
      real    h
C
      real integral     
      real x
      integer i
C
      real f   
C
      integral = (f(local_a) + f(local_b))/2.0
      x = local_a
      do 100 i = 1 , local_n-1 
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
      real function f( x)
      real x
C  Calculate f(x).   
C  Store calculation in return_val.    
      f = x*x
      return
      end
C
