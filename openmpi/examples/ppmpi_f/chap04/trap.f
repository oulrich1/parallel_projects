c  trap.f -- Parallel Trapezoidal Rule, first version
c 
c  Input: None.
c  Output:  Estimate of the integral from a to b of f(x) 
c     using the trapezoidal rule and n trapezoids.
c 
c  Algorithm:
c     1.  Each process calculates "its" interval of 
c         integration.
c     2.  Each process estimates the integral of f(x)
c         over its interval using the trapezoidal rule.
c     3a. Each process != 0 sends its integral to 0.
c     3b. Process 0 sums the calculations received from
c         the individual processes and prints the result.
c 
c  Notes:  
c     1.  f(x), a, b, and n are all hardwired.
c     2.  Assumes number of processes (p) evenly divides
c         number of trapezoids (n = 1024)
c
c  See Chap. 4, pp. 56 & ff. in PPMPI.
c
	program trapezoidal
c
	include 'mpif.h'
c
	integer   my_rank
	integer   p
	real      a
	real      b
	integer   n
	real      h
	real      local_a
	real      local_b
	integer   local_n
	real      integral
	real      total 
	integer   source
	integer   dest
	integer   tag
	integer   status(MPI_STATUS_SIZE)
        integer   ierr
c
        real      Trap
c
	data a, b, n, dest, tag /0.0, 1.0, 1024, 0, 0/
	
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr)
	
	h = (b-a)/n
	local_n = n/p
	
	local_a = a + my_rank*local_n*h
	local_b = local_a + local_n*h
	integral = Trap(local_a, local_b, local_n, h)
	
	if (my_rank .EQ. 0) then
	    total = integral
	    do 100 source = 1, p-1
	        call MPI_RECV(integral, 1, MPI_REAL, source, tag, 
     +              MPI_COMM_WORLD, status, ierr)
	        total = total + integral
 100        continue
	else
	    call MPI_SEND(integral, 1, MPI_REAL, dest, 
     +          tag, MPI_COMM_WORLD, ierr)
        endif
	
	if (my_rank .EQ. 0) then
            write(6,200) n
 200        format(' ','With n = ',I4,' trapezoids, our estimate')
            write(6,300) a, b, total
 300        format(' ','of the integral from ',f6.2,' to ',f6.2,
     +             ' = ',f11.5)
        endif

	call MPI_FINALIZE(ierr) 
        end 
c	
c	
	real function Trap(local_a, local_b, local_n, h)
	real     local_a
	real     local_b
	integer  local_n
	real     h 
c
	real     integral
	real     x 
	real     i 
c
	real     f
c
	integral = (f(local_a) + f(local_b))/2.0 
	x = local_a 
	do 100 i = 1, local_n-1
	   x = x + h 
	   integral = integral + f(x) 
 100    continue
	Trap = integral*h 
        return
        end
c	
c	
	real function f(x)
	real x
	real return_val 
        return_val = x*x
        f = return_val
	return 
	end
