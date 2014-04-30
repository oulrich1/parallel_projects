C  mat_mult.f
C  Multiply a sequence of 2x2 matrices -- 1 factor from 
C      each process.  Erroneous.
C 
C  Input: none
C 
C  Output: product of a sequence of 2x2 matrices
C 
C  Algorithm
C     1. Generate local matrix
C     2. Send local matrix to process 0
C     3  if (my_rank == 0)
C     3b.    for each process, receive matrix and 
C               multiply by product
C     3c.    print product
C 
C  Notes: 
C     1. The matrices are stored as linear arrays.  The 
C        correspondence is row major: Matrix[i][j] <-> 
C        Array[2*i + j]
C     2. Local matrices have the form
C             [my_rank    my_rank+1]
C             [my_rank+2  my_rank  ]
C 
C  See Chap 9, pp. 188 & ff in PPMPI
C 
      PROGRAM MatMult
      INCLUDE 'mpif.h'
      integer  MATRIX_ORDER 
      integer  ARRAY_ORDER
      parameter (MATRIX_ORDER = 2, ARRAY_ORDER = 4)  
      real        my_matrix(0:3)
      real        temp(0:3)
      real        product(0:3)
      data product /1,0,0,1/
C  product is the identity matrix   
      integer     p
      integer     my_rank
      integer     status(MPI_STATUS_SIZE)
      integer     ierr
      integer     i
      character *15 title
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      call Initialize(my_matrix, my_rank)
C
      call MPI_SEND(my_matrix, ARRAY_ORDER, MPI_REAL , 0, 0,
     +     MPI_COMM_WORLD, ierr )
C
      if (my_rank .EQ. 0)  then
           do 100 i = 0, p-1
              call MPI_RECV(temp, ARRAY_ORDER, MPI_REAL ,
     +             MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
     +              status, ierr)
              call Mult(product, temp)
 100      continue
          call Print_matrix('The product is ', product)
      endif
C
      call MPI_FINALIZE(ierr)
      end
C
C
C ********************************************************************  
      subroutine Initialize(my_matrix, my_rank)
      real  my_matrix(0:3)   
      integer    my_rank       
C
      my_matrix(0) = my_rank
      my_matrix(3) = my_rank
      my_matrix(1) = (my_rank + 1)
      my_matrix(2) = (my_rank + 2)
      return
      end  
C
C
C ********************************************************************  
      subroutine Mult(product, factor)
      real  product(0:3)   
      real  factor(0:3)
      integer  MATRIX_ORDER 
      integer  ARRAY_ORDER
      parameter (MATRIX_ORDER = 2, ARRAY_ORDER = 4)  

      integer    i, j, k
      real  temp(0:3)
C
      do 100 i = 0 , MATRIX_ORDER-1
          do 200 j = 0  , MATRIX_ORDER -1
              temp(i * MATRIX_ORDER +j) = 0.0
              do 300 k = 0, MATRIX_ORDER -1
                  temp(i * MATRIX_ORDER +j) =
     +                temp(i * MATRIX_ORDER +j) +
     +                ( product(i *MATRIX_ORDER + k) *
     +                 factor(k*MATRIX_ORDER + j) )
 300          continue
 200     continue
 100  continue
C
      do 400 i = 0, ARRAY_ORDER-1
          product(i) = temp(i)
 400  continue
C
      return
      end
C
C
C ********************************************************************  
      subroutine Print_matrix(title, matrix)
      character *15  title      
      real  matrix(0:3)
      integer  MATRIX_ORDER
      integer  ARRAY_ORDER
      parameter (MATRIX_ORDER = 2, ARRAY_ORDER = 4)
C
      integer i, j
C
      print *, title
      do 100 i = 0 , MATRIX_ORDER-1
          print *, (matrix(i * MATRIX_ORDER +j),
     +                 j = 0, MATRIX_ORDER -1 )
 100  continue
      return
      end
