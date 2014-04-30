C  serial_mat_mult.f -- multiply two square matrices on a single processor
C 
C  Input:
C      n: order of the matrices
C      A,B: factor matrices -- in row-major order
C 
C  Output:
C      C: product matrix
C 
C  See Chap 7, pp. 111 & ff in PPMPI
C 
      PROGRAM SerMlt
      integer       MAX_ORDER
      parameter     (MAX_ORDER = 10)
      integer       n
      real  	    A(10, 10)
      real  	    B(10, 10)
      real  	    C(10, 10) 
C
      print *, 'What''s the order of the matrices?'
      read *, n
C
      call Read_matrix('Enter A', A, n)
      call Print_matrix('A =            ', A, n)
      call Read_matrix('Enter B', B, n)
      call Print_matrix('B =            ', B, n)
      call Serial_matrix_mult(A, B, C, n)
      call Print_matrix('Their product is', C, n)
C
      end
C
C ***************************************************************  
      subroutine Read_matrix(prompt, A, n)
      character  *10  prompt 
      integer   n  
      real      A(n, n)               
      integer i, j
C
      print *,   prompt 
      do 100 i = 1, n 
             read *,  (A(i,j), j = 1,n)
 100  continue
      return
      end  
C
C
C ***************************************************************  
C  MATRIX_T is a two-dimensional array of floats   
      subroutine Serial_matrix_mult(A, B, C, n)
      integer n
      real    A(n,n)    
      real    B(n,n)    
      real    C(n,n)        
C
      integer i, j, k
C
      call Print_matrix('In Ser_mult A = ', A, n)
      call Print_matrix('In Ser_mult B = ', B, n)
C
      do 100 i = 1, n 
          do 200 j =1, n      
              C(i,j) = 0.0
              do 300 k = 1, n  
                  C(i,j) = C(i,j) + A(i,k)*B(k,j)
 300          continue
              print 400 ,'i = ', i ,' j = ', j, 
     +                 '  c_ij = ' , C(i,j)
 400          format (A, I2, A, I2, A, F7.2)
 200      continue
 100  continue
      return         
      end
C
C
C ***************************************************************  
      subroutine Print_matrix(title, C, n)
      integer n
      character *20  title   
      real           C(n,n)        
C       
      integer i, j
C
      print *, title 
      do 100 i = 1,n
          print 300,( C(i, j), j = 1, n)
 300          format (1x, 10F7.2)
 100  continue
      return
            end
