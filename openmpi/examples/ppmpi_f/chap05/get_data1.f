C  get_data1.f -- Parallel Trapezoidal Rule; uses a hand-coded 
C      tree-structured broadcast.
C 
C  Input: 
C     a, b: limits of integration.
C     n: number of trapezoids.
C  Output:  Estimate of the integral from a to b of f(x) 
C     using the trapezoidal rule and n trapezoids.
C 
C  Notes: 
C     1.  f(x) is hardwired.
C     2.  Assumes number of processes (p) evenly divides
C         number of trapezoids (n)
C 
C  See Chap. 5, pp. 65 & ff. in PPMPI.
C
      PROGRAM gd1
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
      data   dest, tag /0, 0/
C
      real Trap
C
C  Let the system do what it needs to start up MPI   
      call 	MPI_INIT(ierr)
C
C  Get my process rank   
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
C  Find out how many processes are being used   
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
C
      call Get_data1( a,  b,  n, my_rank, p)
C
      h = (b-a)/n
      local_n = n/p
C
C Length of each process' interval of 
C integration = local_n*h.  So my interval
C starts at:   
      local_a = a + my_rank*local_n*h
      local_b = local_a + local_n*h
      integral = Trap(local_a, local_b, local_n, h)
C
C  Add up the integrals calculated by each process   
      if (my_rank .EQ. 0)  then
          total = integral
          do 100 source = 1, p-1 
              call MPI_RECV( integral, 1, MPI_REAL , source, tag,
     +             MPI_COMM_WORLD,  status, ierr )
	      total = total + integral
 100      continue
      else  
          call MPI_SEND( integral, 1, MPI_REAL , dest,
     +         tag, MPI_COMM_WORLD, ierr )
      endif
C
C  Print the result   
      if (my_rank .EQ. 0)  then
          print *, 'With n = ', n, ' trapezoids, our estimate'
          print *, ' of the integral from ', a, ' to ', b, ' = ',total
      endif
C
C  Shut down MPI   
      call MPI_FINALIZE(ierr)
      end
C
C
C ******************************************************************  
C Ceiling of log_2(x) is just the number of times
C times x-1 can be divided by 2 until the quotient
C is 0.  Dividing by 2 is the same as right shift.
C 
      integer function Ceiling_log2(x)
      integer x
C  Use unsigned so that right shift will fill
C  leftmost bit with 0
C
      integer temp
      integer result
      temp =  x - 1
      result = 0
C
      do while (temp .NE. 0) 
           temp = temp / 2
           result = result + 1
      end do
      Ceiling_log2 = result
      return 
      end
C
C
C ******************************************************  
      integer function I_receive(stage, my_rank, source)
      integer   stage        
      integer   my_rank      
      integer   source   
      integer   power_2_stage
      integer   i
C
C  2^stage = 1 << stage   
      power_2_stage = 1
      do 100 i = 1,  stage
	   power_2_stage = power_2_stage * 2
 100  continue
      if ( power_2_stage .LE. my_rank .and.
     +     my_rank .LT. 2*power_2_stage) then
          source = my_rank - power_2_stage
          I_receive = 1
      else
	  I_receive = 0
      endif
      return
      end
C
C ***********************************************************  
      integer function I_send(stage, my_rank, p, dest)
      integer   stage      
      integer   my_rank    
      integer   p          
      integer   dest   
      integer   power_2_stage
      integer   i
C
C  2^stage = 1 << stage   
      power_2_stage = 1
      do 100 i = 1, stage
	   power_2_stage = power_2_stage * 2
 100  continue
      if (my_rank .LT. power_2_stage) then
          dest = my_rank + power_2_stage
	  if (dest .GE. p) then
              I_send = 0
	  else
	      I_send = 1
          endif
      else
          I_send = 0
      endif
      return
      end
C
C ******************************************************************  
      subroutine Send(a, b, n, dest)
      real  a      
      real  b      
      integer    n      
      integer    dest
      integer ierr
      include 'mpif.h'  
       
C
      call MPI_SEND( a, 1, MPI_REAL , dest, 0, MPI_COMM_WORLD, ierr )
      call MPI_SEND( b, 1, MPI_REAL , dest, 1, MPI_COMM_WORLD, ierr )
      call MPI_SEND( n, 1, MPI_INTEGER, dest, 2,
     +                   MPI_COMM_WORLD, ierr )
      return
      end
C
C
C ******************************************************************  
      subroutine Receive(a, b, n, source)
      real  a   
      real  b   
      integer    n   
      integer    source
      include 'mpif.h'
C
      integer   status(MPI_STATUS_SIZE)
      integer ierr
     
C
      call MPI_RECV(a, 1, MPI_REAL , source, 0,
     +     MPI_COMM_WORLD,  status, ierr )
      call MPI_RECV(b, 1, MPI_REAL , source, 1,
     +     MPI_COMM_WORLD,  status, ierr )
      call MPI_RECV(n, 1, MPI_INTEGER, source, 2,
     +     MPI_COMM_WORLD,  status, ierr )
      return
      end
C
C
C ******************************************************************  
C Function Get_data1
C Reads in the user input a, b, and n.
C Input parameters:
C     1.  int my_rank:  rank of current process.
C     2.  int p:  number of processes.
C Output parameters:  
C     1.  float* a:  pointer to left endpoint a.
C     2.  float* b:  pointer to right endpoint b.
C     3.  int* n:  pointer to number of trapezoids.
C Algorithm:
C     1.  Process 0 prompts user for input and
C         reads in the values.
C     2.  Process 0 sends input values to other
C         processes using hand-coded tree-structured
C         broadcast.
C 
      subroutine Get_data1(a, b, n, my_rank, p)
      real  a     
      real  b     
      integer    n     
      integer    my_rank   
      integer    p         
C
      integer source
      integer dest
      integer stage
C
      integer Ceiling_log2
      integer I_receive 
      integer I_send 
	  
C
      if (my_rank .EQ. 0) then
          print *, 'Enter a, b, and n'
          read *, a, b, n
      endif
C
      do 100 stage = 0, Ceiling_log2(p)-1  
	  if (I_receive(stage, my_rank, source).EQ. 1) then
              call Receive(a, b, n, source)
          else if (I_send(stage, my_rank, p,  dest).EQ.1)then
	      call Send(a, b, n, dest)
	  endif
 100  continue
      return
      end
C
C
C ******************************************************************  
      real function Trap(local_a, local_b, local_n, h)
      real  local_a    
      real  local_b    
      integer  local_n    
      real  h          
C
      real integral     
      real x
      integer i
C
      real f  
C
      integral = (f(local_a) + f(local_b))/2.0
      x = local_a
      do 100 i = 1 ,local_n-1   
          x = x + h
          integral = integral + f(x)
 100  continue
      Trap = integral*h
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
C
