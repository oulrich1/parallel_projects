C  serial_dot.f -- compute a dot product on a single processor.
C 
C  Input: 
C      n: order of vectors
C      x, y:  the vectors
C 
C  Output:
C      the dot product of x and y.
C 
C  Note:  Arrays containing vectors are statically allocated.
C 
C  See Chap 5, p. 75 in PPMPI.
C 
C
      PROGRAM serdot
      INCLUDE 'mpif.h'
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      real  x(MAX_ORDER)
      real  y(MAX_ORDER)
      integer    n
      real  dot
C
      real Serial_dot
C
      print *, 'Enter the order of the vectors'
      read *, n
      call Read_vector('the first vector', x, n)
      call Read_vector('the second vector', y, n)
      dot = Serial_dot(x, y, n)
      print *, 'The dot product is ', dot
      end
C
C
C ***************************************************************  
      subroutine Read_vector(prompt, v, n)
      integer n
      character *20  prompt   
      real  v(n)              
      integer i
C
      print *, 'Enter ', prompt, 'data (return after each entry): '
      do 100 i = 1, n
          read *,   v(i)
 100  continue
      return
      end
C
C
C ***************************************************************  
      real function Serial_dot(x, y, n)
      integer n
      real  x(n)   
      real  y(n)   
C
      integer    i
      real  sum 
C
      sum = 0
      do 100 i = 1 , n  
          sum = sum + x(i)*y(i)
 100  continue
      Serial_dot = sum
      return 
      end  
