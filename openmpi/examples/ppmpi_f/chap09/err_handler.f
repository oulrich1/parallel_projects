C  err_handler.f -- change default error handler in MPI to 
C      MPI_ERRORS_RETURN
C 
C  Input: none.
C  Output: Error message from each process.
C 
C  See Chap 9, pp. 210 & ff in PPMPI.
C 
      PROGRAM ErrHandler
      INCLUDE 'mpif.h'
      integer p
      integer my_rank
      character *100  error_message
      integer message_length
      integer error_code
      integer x(0:1)
      integer count  
      parameter (count = 2)
      integer    comm
      integer    ierr
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      call MPI_ERRHANDLER_SET(MPI_COMM_WORLD,
     +               MPI_ERRORS_RETURN, ierr)
      error_code = ierr
C
      if (my_rank .EQ. 0)  then
          x(0) = 1
          x(1) = 2
      endif
      call  MPI_BCAST(x, count, MPI_INTEGER, 0, 
     +                   MPI_COMM_WORLD, ierr)
C     call MPI_BCAST(x, count, MPI_INTEGER, 0, comm, ierr)
      error_code = ierr
      if (error_code .NE. MPI_SUCCESS)  then
          call MPI_ERROR_STRING(error_code, error_message,
     +          message_length, ierr)
          print *, 'Error in call to call MPI_BCAST = ',
     +                 error_message
          print *, 'Exiting from function XXX'
          call MPI_ABORT(MPI_COMM_WORLD, -1, ierr )
      endif
C
C
      call MPI_FINALIZE(ierr)
      end
