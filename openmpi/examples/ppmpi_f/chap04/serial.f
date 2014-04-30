C serial.f -- calculate definite integral using trapezoidal rule.
C
C The function f(x) is hardwired.
C Input: a, b, n.
C Output: estimate of integral from a to b of f(x)
C    using n trapezoids.
C 
C See Chapter 4, pp. 53 & ff. in PPMPI.
C
      PROGRAM serial
      INCLUDE 'mpif.h'
      real  integral     
      real  a
      real  b
      integer  n
      real  h            
      real  x
      integer i
C
      real f    
C
      print *, 'Enter a, b, and n'
      read *,  a,  b,  n
C
      h = (b-a)/n
      integral = (f(a) + f(b))/2.0
      x = a
      do 100 i = 1 , n-1
          x = x + h
          integral = integral + f(x)
 100  continue
      integral = integral*h
C
      print *,'With n =', n,' trapezoids, our estimate'
      print *,'of the integral from ', a, ' to ',b, ' = ' , integral
      end
C
C******************************************************
      real function f(x)
      real x
C  Calculate f(x).  Store calculation in return_val.    
C
      f = x*x
      return
      end
