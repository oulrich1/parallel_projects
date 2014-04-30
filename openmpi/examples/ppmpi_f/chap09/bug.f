C  bug.f.  Serial program
C  Program that tries to read a list of floats and sort them in
C     increasing order.
C  Warning!  This program is definitely incorrect!
C 
C  Input:
C     size of list (int)
C     list of floats
C 
C  Output:
C     Sorted list
C 
C  Algorithm:
C     1. Get list size
C     2. Get first input float.
C     3. For each new element
C         (a) read it in
C         (b) use linear search to determine where it
C             should be inserted
C         (c) insert it by shifting greater elements down
C             one
C     4. Print list
C 
C  See Chap 9, pp 180 & ff in PPMPI
C
      PROGRAM Bug
      integer MAX
      parameter (MAX = 100)
      integer num_vals     
      real x(0:99)
      data x /100*0/
      real temp       
      integer i, j, k      
C  values, j:  position to insert     
C  new value.                         
C
      print * ,'How many input values?'
      read *, num_vals
C
      print *, 'Now enter each value: '   
      read *, x(0)
C
      do 100 i = 1, num_vals-1   
          read *, temp
C
C  Determine where to insert   
          j = i - 1
          do while ((j .GT. 0) .AND. (temp .LT. x(j)) )
              j = j - 1
          end do
C
C  Insert
          k = i
C  Warning! Infinite loop when k already > j!!!
          do while ( k .GT. j )
              x(k) = x(k-1)
              k = k + 1
          end do
          x(j) = temp
 100  continue
      end
C
      subroutine Snapshot( title,  num_vals,  x,
     +                       i, j ,  k, temp)
      character *18 title
      integer num_vals
      real    x(0:num_vals-1)
      integer i, j, k
      integer cnt
      real    temp
C
      print *,'*********************'
      print *, title
      print 10, 'num_vals = ',num_vals,', i = ',i,
     +     ', temp = ',temp
 10   format (a, i4, a, i3, a, f10.5)
      print 20, 'j = ',j, ', k = ',k
 20   format (a, i3, a, i3)
      print *, 'x = ', (x(cnt), cnt = 0,i-1)
      print *,'*********************'
      return
      end
C
