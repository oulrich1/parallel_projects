C  cache_test.f -- cache and retrieve a process rank with a communicator
C 
C  Input: none
C  Output: Message from the process whose rank was cached
C 
C  See Chap 8, pp. 139 & ff in PPMPI
C
      PROGRAM Cache
      INCLUDE 'mpif.h'
      integer    p
      integer    my_rank
      integer    io_comm        
      integer    IO_KEY         
      integer    io_rank   
      integer    extra_arg      
      integer    flag
      integer    io_rank_attr
      integer    ierr
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
C  Get a separate communicator for I/O functions by   
C  duplicating MPI_COMM_WORLD                     
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
C
C  Create the attribute key   
      call MPI_KEYVAL_CREATE(MPI_DUP_FN, 
     +      MPI_NULL_DELETE_FN,
     +      IO_KEY, extra_arg,  ierr)
C
C  Set the attribute content   
       io_rank = 0
C
C  Cache the attribute with io_comm   
      call MPI_ATTR_PUT(io_comm, IO_KEY, io_rank, ierr)
C
C  Retrieve the I/O process rank    
      call MPI_ATTR_GET(io_comm, IO_KEY,  io_rank_attr,  
     +                  flag, ierr)
C
C  If flag == 0, something went wrong:   
C  there's no attribute cached.      
      if ((flag .NE. 0) .AND. (my_rank .EQ. io_rank_attr)) then
          print *, 'Greetings from the I/O Process!'
          print *,'My rank is ', io_rank_attr
      endif
C
      call MPI_FINALIZE(ierr)
      end
