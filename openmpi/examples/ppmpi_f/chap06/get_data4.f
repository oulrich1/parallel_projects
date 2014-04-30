C  get_data4.f -- Parallel Trapezoidal Rule.  Uses MPI_Pack/Unpack in
C      distribution of input data.
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
C  See Chap 6., pp. 100 & ff in PPMPI
C 
C
      PROGRAM GetDat
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
      integer   dest    
      integer   tag 
      integer   ierr
      data   dest, tag /0, 0/
C
      real Trap
C
C  Let the system do what it needs to start up MPI   
      call MPI_INIT( ierr )
C
C  Get my process rank   
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
C  Find out how many processes are being used   
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
C
      call Get_data4( a,  b,  n, my_rank)
C
      h = (b-a)/n      
      local_n = n/p    
C
C  Length of each process' interval of 
C  integration = local_n*h.  So my interval
C  starts at:   
      local_a = a + my_rank*local_n*h
      local_b = local_a + local_n*h
      integral = Trap(local_a, local_b, local_n, h)
C
C  Add up the integrals calculated by each process   
      call MPI_REDUCE( integral,  total, 1, MPI_REAL ,
     +     MPI_SUM, 0, MPI_COMM_WORLD, ierr )
C
C  Print the result   
	if (my_rank .EQ. 0) then
            write(6,100) n
 100        format(' ','With n = ',I4,' trapezoids, our estimate')
            write(6,200) a, b, total
 200        format(' ','of the integral from ',f6.2,' to ',f6.2,
     +             ' = ',f11.5)
        endif
C
C  Shut down MPI   
      call MPI_FINALIZE(ierr)
      end
C
C
C ******************************************************************  
      subroutine Get_data4(a, b, n, my_rank)
      real  a     
      real  b     
      integer    n     
      integer    my_rank   
C
      INCLUDE   'mpif.h'
      integer   ierr
      character buffer(100)    
      integer   position       
C      in the buffer             
C
      if (my_rank .EQ. 0) then
          print *,'Enter a, b, and n'
          read *, a, b, n
C
C  Now pack the data into buffer.  Position = 0   
C  says start at beginning of buffer.             
          position = 0
C
C  Position is in/out   
          call MPI_PACK(a, 1, MPI_REAL , buffer, 100,
     +         position, MPI_COMM_WORLD, ierr )
C  Position has been incremented: it now refer-   
C  ences the first free location in buffer.       
C
          call MPI_PACK(b, 1, MPI_REAL , buffer, 100,
     +          position, MPI_COMM_WORLD, ierr )
C  Position has been incremented again.   
C
          call MPI_PACK(n, 1, MPI_INTEGER, buffer, 100,
     +          position, MPI_COMM_WORLD, ierr )
C  Position has been incremented again.   
C
C  Now broadcast contents of buffer   
          call MPI_BCAST(buffer, 100, MPI_PACKED, 0,
     +         MPI_COMM_WORLD, ierr )
       else  
          call MPI_BCAST(buffer, 100, MPI_PACKED, 0,
     +        MPI_COMM_WORLD, ierr )
C
C  Now unpack the contents of buffer   
          position = 0
          call MPI_UNPACK(buffer, 100,  position, a, 1,
     +        MPI_REAL , MPI_COMM_WORLD, ierr )
C  Once again position has been incremented:   
C  it now references the beginning of b.       
C
          call MPI_UNPACK(buffer, 100,  position, b, 1,
     +        MPI_REAL , MPI_COMM_WORLD, ierr )
          call MPI_UNPACK(buffer, 100,  position, n, 1,
     +        MPI_INTEGER, MPI_COMM_WORLD, ierr )
       endif
       return
       end
C
C
C ******************************************************************  
	real function Trap(local_a, local_b, local_n, h)
	real     local_a
	real     local_b
	integer  local_n
	real     h 
C
	real     integral
	real     x 
	real     i 
C
	real     f
C
	integral = (f(local_a) + f(local_b))/2.0 
	x = local_a 
	do 100 i = 1, local_n-1
	   x = x + h 
	   integral = integral + f(x) 
 100    continue
	Trap = integral*h 
        return
        end
C	
C******************************************************************
	real function f(x)
	real x
	real return_val 
        return_val = x*x
        f = return_val
	return 
	end
