C  get_data3.f -- Parallel Trapezoidal Rule.  Builds a derived type
C      for use with the distribution of the input data.
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
C  See Chap 6, pp. 90 & ff in PPMPI
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
      integer    ierr, dest, tag
      data   dest, tag /0, 0/
C
      real Trap     
C
C  Let the system do what it needs to start up MPI   
      call MPI_INIT( ierr )
C
C  Get my process rank   
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr )
C
C  Find out how many processes are being used   
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
C
      call Get_data3( a,  b,  n, my_rank)
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
      subroutine Build_derived_type(a, b, n, 
     +                             mesg_mpi_t)
      INCLUDE       'mpif.h'
      real          a            
      real          b            
      integer       n            
      integer       mesg_mpi_t   
C  pointer to new MPI type   
C
C  The number of elements in each 'block' of the     
C      new type.  For us, 1 each.                    
      integer      block_lengths(3)
      integer      ierr
C  MPI types of the elements.  The 't_i's.'          
      integer       typelist(3)
C
C  Displacement of each element from start of new    
C      type.  The 'd_i's.'                                                   
      integer displacements(3)
C
C  Use for calculating displacements                 
      integer start_address
      integer address
C
      block_lengths(1) = 1
      block_lengths(2) = 1
      block_lengths(3) = 1
C
C  Build a derived datatype consisting of    
C  two floats and an int                     
      typelist(1) = MPI_REAL 
      typelist(2) = MPI_REAL 
      typelist(3) = MPI_INTEGER
C
C  First element, a, is at displacement 0        
      displacements(1) = 0
C
C  Calculate other displacements relative to a   
      call MPI_ADDRESS(a,  start_address, ierr)
C
C  Find address of b and displacement from a     
      call MPI_ADDRESS(b,  address, ierr)
      displacements(2) = address - start_address
C
C  Find address of n and displacement from a     
      call MPI_ADDRESS(n,  address, ierr)
      displacements(3) = address - start_address
C
C  Build the derived datatype   
      call MPI_TYPE_STRUCT(3, block_lengths, displacements,
     +     typelist, mesg_mpi_t, ierr)
C
C  Commit it -- tell system we'll be using it for   
C  communication.                                   
      call MPI_TYPE_COMMIT(mesg_mpi_t, ierr)
      return   
      end  
C
C
C ******************************************************************  
      subroutine Get_data3(a, b, n, my_rank)
      real  a     
      real  b     
      integer     n     
      integer     my_rank   
      integer     mesg_mpi_t   
      integer     ierr
      include     'mpif.h'
C  to 3 floats and an int   
C
      if (my_rank .EQ. 0) then
          print *, 'Enter a, b, and n'
          read *, a, b, n
      endif
C
      call Build_derived_type(a, b, n,  mesg_mpi_t)
      call MPI_BCAST(a, 1, mesg_mpi_t, 0, MPI_COMM_WORLD, ierr )
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
C ******************************************************************  	
	real function f(x)
	real x
	real return_val 
        return_val = x*x
        f = return_val
	return 
	end
