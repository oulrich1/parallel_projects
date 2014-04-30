C  parallel_mat_vect.f -- computes a parallel matrix-vector product.  Matrix
C      is distributed by block rows.  Vectors are distributed by blocks.
C 
C  Input:
C      m, n: order of matrix
C      A, x: the matrix and the vector to be multiplied.  A should be
C          input in row-major order.
C 
C  Output:
C      y: the product vector
C 
C  Notes:  
C      1.  Local storage for A, x, and y is statically allocated.
C      2.  Program treats matrices as one-dimensional arrays in
C          order to override Fortran's column-major storage & make
C          the code more similar to the corresponding C source.
C          The matrix element in row i and column j is stored in
C          array entry (i-1)*No_Of_Cols + j.
C      3.  Assumes number of processes (p) evenly divides both
C          m and n.
C 
C  See Chap 5, p. 78 & ff in PPMPI.
C 
C
      PROGRAM parmat
      INCLUDE 'mpif.h'
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      integer        my_rank
      integer        p
      real           local_A(MAX_ORDER * MAX_ORDER)
      real           global_x(MAX_ORDER)
      real           local_x(MAX_ORDER)
      real           local_y(MAX_ORDER)
      integer        m, n
      integer        local_m, local_n
      integer        ierr
C
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      if (my_rank .EQ. 0)  then
          print *, 'Enter the order of the matrix (m x n)'
          read *,  m,  n
      endif
C 
      call MPI_BCAST( m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
C
      local_m = m/p
      local_n = n/p
C
      call Read_matrix('Enter the matrix', local_A, local_m, n, 
     +                     my_rank, p)
      call Print_matrix('We read', local_A, local_m, n, my_rank, p)
C
      call Read_vector('Enter the vector', local_x, local_n, 
     +                     my_rank, p)
      call Print_vector('We read', local_x, local_n, my_rank, p)
C
      call Parallel_matrix_vector_prod(local_A, m, n, local_x, 
     +     global_x, local_y, local_m, local_n)
      call Print_vector('The product is', local_y, local_m, 
     +     my_rank, p)
C
      call MPI_FINALIZE(ierr)
      end
C
C ********************************************************************  
      subroutine Read_matrix(prompt, local_A, local_m, n, my_rank, p)
      character *20       prompt 
      integer             local_m   
      integer             n          
      real                local_A(local_m * n)          
      integer             my_rank   
      integer             p         
C
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      include 'mpif.h'
      integer ierr
      integer             i, j
      real                temp(MAX_ORDER * MAX_ORDER)
C
C  Fill dummy entries in temp with zeroes   
      do 100 i = 1 , p*local_m 
          do 200 j = n+1 , MAX_ORDER
              temp((i-1)*n + j) = 0.0
 200      continue
 100  continue
C
      if (my_rank .EQ. 0)  then
          print *,  prompt, '(one row at a time):'
          do 300 i = 1, p*local_m 
                  read *, (temp((i-1) *MAX_ORDER + j), j=1,n)
 300      continue
      endif

      call MPI_SCATTER(temp, local_m*MAX_ORDER, MPI_REAL , local_A,
     +     local_m*MAX_ORDER, MPI_REAL , 0, MPI_COMM_WORLD, ierr )
C
      return 
      end   
C
C
C ********************************************************************  
      subroutine Read_vector(prompt, local_x, local_n, my_rank, p)
      character *20  prompt
      integer    local_n       
      real       local_x(local_n)            
      integer    my_rank     
      integer    p           
C
      include 'mpif.h'
      integer ierr
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      integer   i
      real temp(MAX_ORDER)
C
      if (my_rank .EQ. 0)  then
          print *, prompt
          read *, ( temp(i), i = 1, p*local_n)
      endif
      call MPI_SCATTER(temp, local_n, MPI_REAL , local_x, local_n, 
     +           MPI_REAL ,0, MPI_COMM_WORLD, ierr )
C
      return
      end   
C
C
C ********************************************************************  
C  All arrays are allocated in calling program   
      subroutine Parallel_matrix_vector_prod(local_A, m, n, local_x, 
     +                          global_x, local_y, local_m, local_n)
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      integer        m            
      integer        n 
      integer        local_m      
      integer        local_n   
      real           local_A(local_m * MAX_ORDER)                         
      real           local_x(local_n)    
      real           global_x(n)   
      real           local_y(local_m)       
C
C  local_m = m/p, local_n = n/p   
C
      include 'mpif.h'
      integer ierr
      integer i, j
C
      call MPI_ALLGATHER(local_x, local_n, MPI_REAL ,
     +                global_x, local_n, MPI_REAL ,
     +                MPI_COMM_WORLD, ierr )
      do 100 i = 1, local_m 
          local_y(i) = 0.0
          do 200 j = 1, n 
             local_y(i) = 
     +           local_y(i) + local_A((i-1)*MAX_ORDER +j)*global_x(j)
 200      continue
 100  continue
      return 
      end   
C
C
C ********************************************************************  
      subroutine Print_matrix(title, local_A, local_m, n, my_rank, p)
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      character *8        title         
      integer             local_m     
      integer             n
      real                local_A(local_m * MAX_ORDER)           
      integer             my_rank     
      integer             p           
C
      include 'mpif.h'
      integer ierr
      integer   i, j
      real temp(MAX_ORDER * MAX_ORDER)
C
      call MPI_GATHER(local_A, local_m*MAX_ORDER, MPI_REAL , temp,
     +     local_m*MAX_ORDER, MPI_REAL , 0, MPI_COMM_WORLD, ierr )
C
      if (my_rank .EQ. 0)  then
          print *, title
          do 100 i = 1, p*local_m
                  print *,  (temp((i-1)*MAX_ORDER + j), j=1, n)
 100      continue
      endif
C	      
      return
      end
C
C
C ********************************************************************  
      subroutine Print_vector(title, local_y, local_m, my_rank, p)
      character *8  title 
      integer    local_m
      real       local_y(local_m)        
      integer    my_rank     
      integer    p           
C
      include 'mpif.h'
      integer ierr
      integer   i
      integer MAX_ORDER
      parameter (MAX_ORDER = 100)
      real temp(MAX_ORDER)
C
      call MPI_GATHER(local_y, local_m, MPI_REAL, temp, local_m, 
     +     MPI_REAL,0 , MPI_COMM_WORLD, ierr )
C
      if (my_rank .EQ. 0)  then
          print *,  title
          do 100 i = 1, local_m*p
              print *, temp(i)
 100      continue
          print *,' ' 
      endif
C
      return
      end   
