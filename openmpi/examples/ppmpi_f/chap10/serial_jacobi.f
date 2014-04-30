C  serial_jacobi.f -- serial version of Jacobi's method for solving
C      the linear system Ax = b.
C 
C  Input:
C      n:  order of system
C      tol:  convergence tolerance
C      max_iter:  maximum number of iterations
C      A:  coefficient matrix
C      b:  right-hand side of system
C 
C  Output:
C      x:  the solution if the method converges
C      max_iter:  if the method fails to converge
C 
C  Notes:  
C      1.  A should be strongly diagonally dominant in
C          order to insure convergence.
C 
C  See Chap 10, pp. 218 & ff in PPMPI.
C 
C
      PROGRAM SerJacobi
      real     A(0:11, 0:11)
      real     x(0:11)
      real     b(0:11)
      integer  MAX_DIM
      parameter (MAX_DIM = 12)
      integer  n
      real     tol
      integer  max_iter
      integer  converged
      integer  Jacobi
C
      print *,'Enter n, tolerance, and max number 
     + of iterations'
      read *,  n,  tol,  max_iter 
      call Read_matrix('Enter the matrix          ', A, n)
      call Read_vector('Enter the right-hand side ', b, n)
C
      converged = Jacobi(A, x, b, n, tol, max_iter)
C
      if (converged .EQ. 1 ) then
          call Print_vector('The solution is           ', x, n)
      else
          print *, 'Failed to converge in ', max_iter,
     +                   ' iterations'
      endif
      end
C
C
C *******************************************************************  
C  Return 1 if iteration converged, 0 otherwise   
C  MATRIX_T is just a 2-dimensional array         
      integer function Jacobi(A, x, b, n, tol, max_iter)
      real     A(0:n-1, 0:n-1)          
      real     x(0:n-1)        
      real     b(0:n-1)        
      integer  n          
      real     tol        
      integer    max_iter   
      integer    i, j
      integer    iter_num
      real  retval
      real  x_old(0:11)
C
      real Distance 
C
C  Initialize x   
      do 100 i = 0 ,n-1
          x(i) = b(i)
 100  continue
C
      iter_num = 0
      retval = 99 
      do while ((iter_num .LT. max_iter) .AND.
     +          (retval .GE. tol))
          iter_num = iter_num + 1
C
          do 200 i = 0 ,n-1
              x_old(i) = x(i)
 200      continue
C
          do 300 i = 0, n-1
              x(i) = b(i)
              do 400 j = 0 , i-1
                  x(i) = x(i) - A(i, j)*x_old(j)
 400          continue
              do 500 j = i+1, n-1
                  x(i) = x(i) - A(i, j)*x_old(j)
 500          continue
              x(i) = x(i)/A(i,i)
 300      continue
          retval = Distance(x,x_old,n)
      end do
C
      if (retval .LT. tol) then
          Jacobi = 1
      else
          Jacobi = 0
      endif
C
      return
      end 
C
C
C *******************************************************************  
      real function Distance(x, y, n)
      real x(0:n-1)
      real y(0:n-1)
      integer n
      integer i
      real sum
C
      sum = 0.0
      do 100 i = 0 , n - 1
          sum = sum + ((x(i) - y(i))*(x(i) - y(i)))
 100  continue
      Distance =  sqrt(sum)
      return
      end
C
C
C *******************************************************************  
      subroutine Read_matrix(prompt, A, n)
      character *27     prompt   
      real    A(0:n-1, 0:n-1)        
      integer n        
      integer i, j
C
      print *, prompt 
      do 100 i = 0  ,n-1
          read *,( A(i,j), j = 0, n-1)
 100  continue
      return
      end
C
C
C *******************************************************************  
      subroutine Read_vector(prompt, x, n)
      character *27  prompt   
      real  x(0:n-1)      
      integer    n        
      integer i
C
      print *, prompt
      read *, (x(i), i = 0, n-1)
      return
      end   
C
C
C *******************************************************************  
      subroutine Print_matrix(title, A, n)
      character *27    title
      real    A(0:n-1, 0:n-1)        
      integer n        
      integer i, j
C
      print *, title
      do 100 i = 0  ,n-1
          print *,( A(i,j), j = 0, n-1)
 100  continue
      return
      end
C
C***********************************************************************
      subroutine Print_vector(title, x, n)
      character *27  title   
      real  x(0:n-1)      
      integer    n        
      integer i
C
      print *, title
      print *, (x(i), i = 0, n-1)
      return
      end   
