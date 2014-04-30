C  serial_bitonic.f -- serial bitonic sort of randomly generated list
C      of integers
C 
C  Input:
C      n: the length of the list -- must be a power of 2.
C 
C  Output:
C      The unsorted list and the sorted list.
C 
C  Notes:
C      1.  The list is statically allocated -- size specified in MAX.
C      2.  Keys are in the range 0 -- KEY_MAX-1.
C 
C  See Chap 14, pp. 316 & ff. in PPMPI.
C
C *******************************************************************  
      PROGRAM SerBit
      integer    list_length
      integer    n
      integer    start_index
      integer    ordering
      integer    A(0:16383)
      integer    KEY_MAX
      common     KEY_MAX
C  Successive subsequences will switch between
C  increasing and decreasing bitonic splits.
      integer    INCR
      integer    DECR
      parameter  (INCR=0,  DECR= 1)
C   
      integer    Reverse
C
      KEY_MAX = 16384
C
      print *,'Enter the list size (a power of 2)'
      read *,  n
C
      call Generate_list(n, A)
C 
      call Print_list('The unsorted list is', n, A)
C 
      list_length = 2
      do while (list_length .LE. n)
          ordering = INCR
          start_index = 0
          do while (start_index .LT. n )  
              if (ordering .EQ. INCR) then
                  call Bitonic_sort_incr(list_length,
     +                  A(start_index) )
              else
                  call Bitonic_sort_decr(list_length,
     +                  A(start_index) )
              endif
              start_index = start_index + list_length
              ordering = Reverse(ordering)
          end do
          list_length = list_length*2
      end do
C
      call Print_list('The sorted list is ', n, A)
C
      end
C
C
C *******************************************************************  
      subroutine Generate_list( n,  A) 
      integer n
      integer A(n)
C
      integer KEY_MAX
      common KEY_MAX
      integer i
      integer iseed
      integer Random
C
      iseed = 20029
C
      do 100 i = 1  , n   
          Random = (i * iseed)        
          A(i) = MOD (Random, KEY_MAX)
 100  continue
      return
      end   
C
C
C *******************************************************************  
      subroutine Bitonic_sort_incr(length, B)
      integer     length   
      integer     B(0:length-1)      
C
      integer i
      integer half_way
C
C  This is the bitonic split   
      half_way = length/2
      do 100 i = 0 , half_way-1 
          if (B(i) .GT. B(half_way + i)) then
              call Swap( B(i), B(half_way+i ) )
          endif
 100  continue
C
      if (length .GT. 2)  then
          call Bitonic_sort_incr(length/2, B)
          call Bitonic_sort_incr(length/2, B(half_way) )
      endif
C
      return
      end
C
C
C *******************************************************************  
      subroutine Bitonic_sort_decr(length, B)
      integer    length    
      integer    B(0:length-1)
C       
      integer i
      integer half_way
C
C  This is the bitonic split   
      half_way = length/2
      do 100 i = 0 , half_way-1  
          if (B(i) .LT. B(half_way + i)) then
              call Swap( B(i), B(half_way+i) )
          endif
 100  continue
C
      if (length .GT. 2)  then
          call Bitonic_sort_decr(length/2, B)
          call Bitonic_sort_decr(length/2, B(half_way) )
      endif
C
      return
      end
C
C *******************************************************************  
        subroutine Print_list(title,  n, A) 
        character *20 title
        integer    n
        integer    A(n)
C
        integer i
C
        print *, title
        do 100 i = 1 , n 
          print 200,  A(i)
 200      format (I7)
 100    continue
        return
        end 
C *******************************************************************   
        integer function Reverse(ordering)
C
        integer ordering
C
C  Successive subsequences will switch between
C  increasing and decreasing bitonic splits.
        integer    INCR
        integer    DECR
        parameter  (INCR = 0,  DECR = 1)
C
        if ( ordering .EQ. INCR) then
            Reverse =  DECR
        else
            Reverse =  INCR
        endif
        return 
        end
C *******************************************************************   
      subroutine Swap(a, b)
C
      integer a
      integer b
      integer temp
C
      temp = a
      a = b
      b = temp
C
      return 
               end
