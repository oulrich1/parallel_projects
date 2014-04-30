C  ciof.h -- header file for use with cio.f -- basic collective
C      I/O functions
C
C  See Chap 8, pp. 142 & ff in PPMPI
C
C  Functions:
      integer Cache_io_rank
      integer Copy_attr
      integer Get_io_rank
      integer Cread
      integer Cprint 
      integer Cerror_test
C
C  All process ranks < HUGE   
      integer  HUGE
      integer  NO_IO_ATTR
      parameter (HUGE=32768, NO_IO_ATTR= -1) 
C
C  with C progs BUFSIZ is defined in stdio.h   
      integer   BUFSIZE
      parameter  (BUFSIZE = 1000)
      character  *1000  io_buf
      common  /CioComm/   io_buf
C
      integer   error_buf
      integer   error_bufsiz
      parameter (error_bufsiz = 0)
C
C  Key identifying IO_Attribute   
      integer    IO_KEY
      common  /CIOKey/   IO_KEY
C  Be sure to add the following line to the main
C  program which uses any CIO functions!
C       IO_KEY = MPI_KEYVAL_INVALID
C 
      logical DEBUG
      parameter (DEBUG = .FALSE.)  
C
      integer extra_arg 
      common /EArg/ extra_arg
C  End of io.h   
