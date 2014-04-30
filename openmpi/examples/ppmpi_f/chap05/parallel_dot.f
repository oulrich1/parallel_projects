C  parallel_dot.f -- compute a dot product of a vector distributed among
C      the processes.  Uses a block distribution of the vectors.
C 
C  Input: 
C      n: global order of vectors
C      x, y:  the vectors
C 
C  Output:
C      the dot product of x and y.
C 
C  Note:  Arrays containing vectors are statically allocated.  Assumes
C      n, the global order of the vectors, is divisible by p, the number
C      of processes.
C 
C  See Chap 5, pp. 75 & ff in PPMPI.
C 
      PROGRAM ParDot
      INCLUDE 'mpif.h'
      integer    MAX_LOCAL_ORDER
      parameter  (MAX_LOCAL_ORDER = 100)
      real       local_x(MAX_LOCAL_ORDER)
      real       local_y(MAX_LOCAL_ORDER)
      integer    n
      integer    n_bar    
      real       dot
      integer    p
      integer    my_rank
      integer ierr
C   
      real Parallel_dot
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      if (my_rank .EQ. 0)  then
          print *,'Enter the order of the vectors: '
          read *, n
      endif
      call MPI_BCAST(n, 1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      n_bar = n/p
C
      call Read_vector('the first vector',local_x, n_bar, p, my_rank)
      call Read_vector('the second vector',local_y, n_bar, p,my_rank)
C
      dot = Parallel_dot(local_x, local_y, n_bar)
C
      if (my_rank .EQ. 0) then 
	  print *, 'The dot product is: ', dot
      endif
C
      call MPI_FINALIZE(ierr)
      end
C
C
C ***************************************************************  
      subroutine Read_vector(prompt, local_v, n_bar, p, my_rank)
      include 'mpif.h'
      integer    MAX_LOCAL_ORDER
      parameter  (MAX_LOCAL_ORDER = 100)
      character *20  prompt
      integer    n_bar      
      real       local_v(n_bar)          
      integer    p           
      integer    my_rank     
      integer    i, q
      real      temp(MAX_LOCAL_ORDER)
      integer   status(MPI_STATUS_SIZE)
      integer   ierr     
C
      if (my_rank .EQ. 0)  then
          print *, 'Enter ', prompt, ' data(return after each entry)'
          do 100 i = 1 , n_bar
              read *,  local_v(i)
 100      continue
          do 200 q = 1 , p-1  
              do 300 i = 1, n_bar 
                  read *, temp(i)
 300          continue
              call MPI_SEND(temp(1), n_bar, MPI_REAL , q, 0, 
     +                    MPI_COMM_WORLD, ierr )
 200      continue
      else  
          call MPI_RECV(local_v(1), n_bar, MPI_REAL , 0, 0, 
     +                    MPI_COMM_WORLD, status, ierr)
      endif
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
      sum = 0
C
      do 100 i = 1, n 
          sum = sum + x(i)*y(i)
 100  continue
      Serial_dot = sum
      return 
      end  
C
C ***************************************************************  
      real function Parallel_dot(local_x, local_y, n_bar)
      include 'mpif.h'
      integer n_bar
      real  local_x(n_bar)   
      real  local_y(n_bar)   
C
      real  local_dot
      real  dot
      integer ierr
C  
      real  Serial_dot
      dot = 0
C
      local_dot = Serial_dot(local_x, local_y, n_bar)
      call MPI_REDUCE( local_dot,  dot, 1, MPI_REAL ,
     +                  MPI_SUM, 0, MPI_COMM_WORLD, ierr )
      Parallel_dot = dot
      return 
       end
