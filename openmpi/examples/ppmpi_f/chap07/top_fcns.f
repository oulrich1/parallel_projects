C  top_fcns.f -- test basic topology functions
C 
C  Input: none
C  Output: results of calls to various functions testing topology
C      creation
C 
C  Algorithm:
C      1.  Build a 2-dimensional Cartesian communicator from
C          MPI_Comm_world
C      2.  Print topology information for each process
C      3.  Use MPI_Cart_sub to build a communicator for each
C          row of the Cartesian communicator
C      4.  Carry out a broadcast across each row communicator
C      5.  Print results of broadcast
C      6.  Use MPI_Cart_sub to build a communicator for each
C          column of the Cartesian communicator
C      7.  Carry out a broadcast across each column communicator
C      8.  Print results of broadcast
C 
C  Note: Assumes the number of process, p, is a perfect square
C 
C  See Chap 7, pp. 121 & ff in PPMPI
C 
      PROGRAM TopFcn
      INCLUDE 'mpif.h'
      integer       p
      real          preal
      integer       my_rank
      integer       q
      integer       grid_comm
      integer       dim_sizes(0:1)
      integer       wrap_around(0:1)
      integer       reorder
      integer       coordinates(0:1)
      integer       my_grid_rank
      integer       grid_rank
      integer       free_coords(0:1)
      integer       row_comm
      integer       col_comm
      integer       row_test
      integer       col_test
      integer       ierr
C
C
      reorder = 1
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      preal = p
      q = sqrt(preal)
C
      dim_sizes(0) = q
      dim_sizes(1) = q
      wrap_around(0) = 1
      wrap_around(1) = 1
      call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dim_sizes,
     +      wrap_around, reorder,  grid_comm, ierr)
C
      call MPI_COMM_RANK(grid_comm,  my_grid_rank, ierr)
      call MPI_CART_COORDS(grid_comm, my_grid_rank, 2,
     +         coordinates, ierr)
C
      call MPI_CART_RANK(grid_comm, coordinates, grid_rank,
     +         ierr)
C
      print 100,'Proc' ,my_rank, '> grid_rank=',my_grid_rank,
     +        ', coords = (', coordinates(0), ',', coordinates(1),
     +        ') grid_rank = ',grid_rank
 100  format (a, i3, a, i3, a, i2, a, i2, a, i3)
C
      free_coords(0) = 0
      free_coords(1) = 1
      call MPI_CART_SUB(grid_comm, free_coords,  row_comm, ierr)
      if (coordinates(1) .EQ. 0)then
          row_test = coordinates(0)
      else
          row_test = -1
      endif
      call MPI_BCAST(row_test, 1,MPI_INTEGER, 0,row_comm, ierr)
      print 200 ,'Proc', my_rank, ' > coords = (',
     +         coordinates(0), ',', coordinates(1),
     +        ') row_test = ', row_test
 200  format (a, i3, a, i2, a, i2, a,i3)     
C
      free_coords(0) = 1
      free_coords(1) = 0
      call MPI_CART_SUB(grid_comm, free_coords,  col_comm, ierr)
      if (coordinates(0) .EQ. 0)then
          col_test = coordinates(1)
      else
          col_test = -1
      endif
      call MPI_BCAST( col_test, 1, MPI_INTEGER,0, col_comm,ierr)
      print 300,'Proc', my_rank, ' > coords = (',
     +         coordinates(0), ',', coordinates(1),
     +        ') col_test = ', col_test
 300  format (a, i3, a, i2, a, i2, a,i3) 
C 
      call MPI_FINALIZE(ierr)
      end
