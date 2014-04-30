C  stdin_test.f -- test whether an MPI implementation allows input from
C      stdin
C 
C  Input: an int
C  Output: prompt for input and int read (if stdin is OK)
C 
C  See Chap 8, p. 154 in PPMPI.
C 
      PROGRAM Access
      INCLUDE 'mpif.h'
      integer x
      integer ierr
C
      call MPI_INIT( ierr)
C
      print *, 'Enter an integer'
      read *,  x
      print *, 'We read x = ', x
C
      call MPI_FINALIZE(ierr)
      end
