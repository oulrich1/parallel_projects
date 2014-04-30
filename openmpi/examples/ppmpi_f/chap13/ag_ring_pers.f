C  ag_ring_pers.f -- ring allgather using persistent communication requests
C 
C  Input: series of blocksizes for allgather, 0 to stop.
C  Output: Contents of gathered array on each process -- list of
C      process ranks, each rank appearing in a block of size blocksize.
C 
C  Note:  array sizes are hardwired in MAX, LOCAL_MAX, and MAX_BYTES.
C 
C  See Chap 13, pp. 301 & ff, in PPMPI.
C
C
C ******************************************************************  
      PROGRAM RingPers
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
      integer       MAX 
      integer       LOCAL_MAX 
      parameter     (MAX = 1024, LOCAL_MAX = 1024)
      integer       p
      integer       my_rank
      real          x(0:1023)
      real          y(0:1023)
      integer       blocksize
      integer       io_comm
      integer       i
      integer       rtnval
      integer       ierr
      character *40 prompt,title
C
      IO_KEY = MPI_KEYVAL_INVALID    
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_DUP(MPI_COMM_WORLD, io_comm, ierr )
      rtnval = Cache_io_rank(MPI_COMM_WORLD, io_comm  )
C
      prompt = 'Enter the local array size (0 to quit)'
      rtnval = Cread(io_comm, 1, prompt, blocksize)
C
      do while(blocksize .GT. 0) 
          do 100 i = 0  , blocksize-1
              x(i) =  my_rank
 100      continue
          call Allgather(x, blocksize, y, MPI_COMM_WORLD  )
          call Print_arrays(io_comm,  'GATHERed_arrays ', 
     +                      y, blocksize)
C  Enter 0 to stop.
          prompt =  'Enter the local array size (0 to quit)'   
          rtnval = Cread(io_comm, 1, prompt, blocksize) 
      end do
C
      call MPI_FINALIZE(ierr)
      end
C
C
C ******************************************************************  
      subroutine Print_arrays(io_comm, title, y, blocksize)
      integer         io_comm     
      character *17     title       
      real            y(128)         
      integer         blocksize   
C
      include     'mpif.h'
      include     'ciof.h'
      integer     ierr
      integer  i, j
      integer  p
      integer  retval
      integer  iitem
      integer  list(128)
      data     list /128*0/
C
      call MPI_COMM_SIZE(io_comm,  p, ierr)
      
      do 100 i = 1 , blocksize*p
          iitem = int(y(i))
          list(i) = iitem
 100  continue
      retval = Cprint(io_comm, blocksize*p, title,list)
      return
      end
C
C
C ******************************************************************  
      subroutine Allgather(x, blocksize, y, ring_comm)
      real      x(0:blocksize)         
      integer   blocksize   
      real      y(0:1023)         
      integer   ring_comm   
C
      INCLUDE    'mpif.h'
      integer    sizeof_real
      parameter  (sizeof_real = 16)
      integer    MAX_BYTES
      parameter  (MAX_BYTES = 1024*sizeof_real)
      integer    i, p, my_rank
      integer    successor, predecessor
      integer    send_offset, recv_offset
      integer    status(MPI_STATUS_SIZE)
      integer    ierr
      integer    send_request
      integer    recv_request
      real       send_buf(0:1023)
      real       recv_buf(0:1023)
      integer    position
C
      call MPI_COMM_SIZE(ring_comm,  p, ierr)
      call MPI_COMM_RANK(ring_comm,  my_rank, ierr)
C
C  Copy x into correct location in y   
      do 100 i = 0, blocksize -1
          y(i + my_rank*blocksize) = x(i)
 100   continue
C
      successor = MOD((my_rank + 1) , p)
      predecessor = MOD((my_rank - 1 + p) , p)
C
      call MPI_SEND_INIT(send_buf, 
     +     blocksize*sizeof_real,
     +     MPI_PACKED, successor, 0, ring_comm,
     +      send_request, ierr)
      call MPI_RECV_INIT(recv_buf, 
     +     blocksize*sizeof_real,
     +     MPI_PACKED, predecessor, 0, ring_comm,
     +     recv_request, ierr )
C
      send_offset = my_rank*blocksize
      do 200 i = 0 , p - 2   
          position = 0
          call MPI_PACK(y(send_offset), blocksize, 
     +         MPI_REAL , send_buf,
     +         MAX_BYTES,  position, ring_comm, ierr)
          call MPI_START( send_request, ierr)
          call MPI_START( recv_request, ierr)
          recv_offset = 
     +        MOD((my_rank - i - 1 + p) , p)*blocksize
          send_offset = recv_offset
          position = 0
          call MPI_WAIT( send_request,  status, ierr)
          call MPI_WAIT( recv_request,  status, ierr)
          call MPI_UNPACK(recv_buf, MAX_BYTES,  
     +        position,
     +        y(recv_offset), blocksize, 
     +        MPI_REAL , ring_comm, ierr)
 200  continue
      return
      end
C
