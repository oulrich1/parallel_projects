C  parallel_jacobi.f -- parallel implementation of Jacobi's method
C      for solving the linear system Ax = b.  Uses block distribution
C      of vectors and block-row distribution of A.
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
C      2.  A is stored in row-major order as a one-dimensional
C          array in order to facilitate block-row distribution
C      3.  In order to make the conversion from C to Fortran
C          easier, array subscripts start at 0 rather than 1
C 
C  See Chap 10, pp. 220 & ff in PPMPI.
C
      PROGRAM ParJacobi
      INCLUDE 'mpif.h'
      integer   p
      integer   my_rank, ierr
      real      A_local(0:12*12-1)
      real      x_local(0:11)
      real      b_local(0:11)
      integer   n
      real      tol
      integer   max_iter
      integer   converged
      integer  MAX_DIM
      parameter (MAX_DIM = 12)
      integer   Parallel_jacobi
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      if (my_rank .EQ. 0)  then
            print *,'Enter n, tolerance, and max number 
     + of iterations'
            read *,  n,  tol,  max_iter
      endif 
      call MPI_BCAST( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( tol, 1, MPI_REAL , 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( max_iter, 1, MPI_INTEGER, 0,
     +                   MPI_COMM_WORLD, ierr )
C
      call Read_matrix('Enter the matrix           ',
     +                          A_local, n, my_rank, p)
      call Read_vector('Enter the right-hand side  ',
     +                         b_local, n, my_rank, p)
C
      converged = Parallel_jacobi(A_local, x_local, b_local, n,
     +     tol, max_iter, p, my_rank)
C
      if (converged .EQ. 1) then
          call Print_vector('The solution is           ',
     +                     x_local, n, my_rank, p)
      else
          if (my_rank .EQ. 0) then
              print *, 'Failed to converge in ', 
     +                  max_iter,' iterations'
          endif
      endif
C
      call MPI_FINALIZE(ierr)
      end
C *******************************************************************  
C  Return 1 if iteration converged, 0 otherwise   
C  MATRIX_T is a 2-dimensional array              
      integer function  Parallel_jacobi(A_local, x_local,
     +           b_local, n, tol, max_iter, p, my_rank)
      real     A_local(0:n*12-1)     
      real     x_local(0:n-1)   
      real     b_local(0:n-1)   
      integer    n           
      real       tol         
      integer    max_iter    
      integer    p           
      integer    my_rank
      INCLUDE 'mpif.h'
      integer    i_local, i_global, j
      integer    n_bar
      integer    iter_num
      real   x_temp1(0:11)
      data   x_temp1 /12*0.0/
      real   x_old(0:11)
      real   x_new(0:11)
      real   retval
      integer  i,  ierr
      integer  MAX_DIM
      parameter (MAX_DIM = 12)
C function 
      real Distance
C
      n_bar = n/p
C
C  Initialize x   
      call MPI_ALLGATHER(b_local, n_bar, MPI_REAL , x_temp1,
     +     n_bar, MPI_REAL , MPI_COMM_WORLD, ierr )
      do 50 i = 0, 11
          x_new(i) = x_temp1(i)
 50   continue
C
      retval = 99
      iter_num = 0
      do while ((iter_num .LT. max_iter) .AND.
     +          (retval .GE. tol))
          iter_num = iter_num + 1
C
C  Interchange x_old and x_new   
          call Swap(x_old, x_new, n)
          do 100 i_local = 0 , n_bar -1
              i_global = i_local + my_rank*n_bar
              x_local(i_local) = b_local(i_local)
              do 200 j = 0  ,i_global-1
                  x_local(i_local) = x_local(i_local) -
     +                 A_local(i_local*MAX_DIM +j)*x_old(j)
 200          continue
              do 300 j = i_global+1,n-1
                  x_local(i_local) = x_local(i_local) -
     +                 A_local(i_local*MAX_DIM + j)*x_old(j)
 300          continue
              x_local(i_local) = x_local(i_local)/
     +                A_local(i_local*MAX_DIM +i_global)
 100      continue
C
          call MPI_ALLGATHER(x_local,n_bar,MPI_REAL,x_new,
     +         n_bar, MPI_REAL , MPI_COMM_WORLD, ierr )
          retval = Distance(x_new,x_old,n)
      end do
C
      if (retval .LT. tol) then
          Parallel_jacobi= 1
      else
          Parallel_jacobi= 0
      endif
C
      return
      end
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
C *******************************************************************  
      subroutine Read_matrix(prompt, A_local, n, my_rank, p)
      character *27    prompt    
      real    A_local(0:n*n-1)   
      integer       n         
      integer       my_rank   
      integer       p         
C
      INCLUDE 'mpif.h'
      integer       i, j, ierr
      real          temp(0:12*12-1)
      integer       n_bar, MAX_DIM
      parameter (MAX_DIM = 12)
C
      n_bar = n/p
C
      do 100 i = 0, 11
         do 200 j = 0, 11
              temp(i*MAX_DIM +j) = 0.0
 200       continue
 100  continue
C
      if (my_rank .EQ. 0)  then
          print *, prompt
          do 300 i = 0, n-1
              read * , (temp(i*MAX_DIM +j), j=0, n-1)
 300      continue
      endif
      call MPI_SCATTER(temp,n_bar*MAX_DIM, MPI_REAL , 
     +    A_local, n_bar*MAX_DIM, MPI_REAL , 0, 
     +    MPI_COMM_WORLD, ierr )
C
      return
      end   
C
C *******************************************************************  
      subroutine Read_vector(prompt, x_local, n, my_rank, p)
      character *27  prompt      
      real  x_local(0:n-1)   
      integer    n           
      integer    my_rank     
      integer    p           
C
      INCLUDE 'mpif.h'
      integer   i, ierr
      real temp(0:11)
      data temp /12*0.0/
      integer   n_bar
C
      n_bar = n/p
C
      if (my_rank .EQ. 0)  then
          print *, prompt
          read *, (temp(i), i = 0, n-1)
      endif
      call MPI_SCATTER(temp, n_bar, MPI_REAL , 
     +        x_local, n_bar, MPI_REAL ,
     +        0, MPI_COMM_WORLD, ierr )
C
      return
      end   
C
C *******************************************************************  
      subroutine Print_matrix(title, A_local, n, my_rank, p)
      character *27    title   
      real    A_local(0:n*n-1)   
      integer       n         
      integer       my_rank   
      integer       p         
C
      INCLUDE 'mpif.h'
      integer       i, j  , ierr
      real          temp(0:12*12-1)
      integer       n_bar, MAX_DIM
      parameter (MAX_DIM = 12)
C
      n_bar = n/p
C
      call MPI_GATHER(A_local, n_bar*MAX_DIM, 
     +        MPI_REAL , temp,  n_bar*MAX_DIM,
     +        MPI_REAL , 0, MPI_COMM_WORLD, ierr )
C
      if (my_rank .EQ. 0)  then
          print *, title
          do 300 i = 0, n-1
              print * , (temp(i*n+j), j=0, n-1)
 300      continue
      endif
C
      return
      end   

C
C *******************************************************************  
      subroutine Print_vector(title, x_local, n, my_rank, p)
      character *27  title   
      real  x_local(0:n-1)   
      integer    n           
      integer    my_rank     
      integer    p           
C
      INCLUDE 'mpif.h'
      integer   i, ierr
      real temp(0:11)
      integer   n_bar
C
      n_bar = n/p
C
      call MPI_GATHER(x_local, n_bar, MPI_REAL , 
     +    temp, n_bar, MPI_REAL ,
     +     0, MPI_COMM_WORLD, ierr )
 
      if (my_rank .EQ. 0)  then
          print *, title
          print *, (temp(i), i = 0, n-1)
      endif
C
      return
      end   
C
C*****************************************************************
      subroutine Swap(a, b, n)
      real a(n), b(n)
      integer n
      real temp(12)
      integer i
      do 100 i = 1, n
          temp(i) = a(i)
          a(i) = b(i)
          b(i) = temp(i)
 100  continue
C
      return
      end
