C  serial_mat_vect.f -- computes a matrix-vector product on a single processor.
C 
C  Input:
C      m, n: order of matrix
C      A, x:  the matrix and the vector to be multiplied
C 
C  Output:
C      y: the product vector
C 
C  Note:  A, x, and y are statically allocated.
C 
C  See Chap 5, p. 78 & ff in PPMPI.
C 
      PROGRAM sermat
      INCLUDE 'mpif.h'
      integer  MAX_ORDER
      parameter (MAX_ORDER = 100)
      real     A(MAX_ORDER, MAX_ORDER)
      real     x(MAX_ORDER)
      real     y(MAX_ORDER)
      integer       m, n
C
      print *,"Enter the order of the matrix (m x n)"
      read *,  m,  n
      call Read_matrix("the matrix", A, m, n)
      call Read_vector("the vector", x, n)
      call Serial_matrix_vector_prod(A, m, n, x, y)
      call Print_vector(y, m)
      end
C
C
C ***************************************************************  
      subroutine Read_matrix(prompt, A, m, n)
      integer   m, n
      character *12  prompt   
      real      A(m, n)            
      integer i, j
C
      print *, "Enter ", prompt, " data(one row at a time): "
      do 100 i = 1, m
              read *, ( A(i,j), j = 1, n)
 100  continue
      return 
      end
C
C
C ***************************************************************  
      subroutine Read_vector(prompt, v, n)    
      character *12  prompt   
      integer    n 
      real       v(n)  
      integer    i
C
      print * , "Enter ", prompt, " data(return after each entry): "
      read *,   (v(i), i = 1,n)
 
C
      return
      end
C
C ***************************************************************  
      subroutine Serial_matrix_vector_prod(A, m, n, x, y)
      integer       m     
      integer       n     
      real     A(m, n)
      real     x(n)   
      real     y(m)   
C
      integer k, j
C
      do 100 k = 1, m 
          y(k) = 0.0
          do 200  j = 1, n
              y(k) = y(k) + A(k,j)*x(j)
 200      continue
 100  continue
C
      return
      end
C
C ***************************************************************  
      subroutine Print_vector(y, n)
      integer n
      real  y(n)   
C   
      integer i
C
      print * , "Result is "
      do 100 i = 1, n
          print *,  y(i)
 100  continue
      print *, "  "
      return
      end   
