C  ag_cube_nblk.f -- hypercube allgather using nonblocking sends and receives
C 
C  Input: series of blocksizes for allgather, 0 to stop.
C  Output: Contents of gathered array on each process -- list of
C      process ranks, each rank appearing in a block of size blocksize.
C 
C  Note:  array sizes are hardwired in MAX and LOCAL_MAX.
C 
C  See Chap 13, pp. 299 & ff, in PPMPI.
C 
C
C ******************************************************************  
      PROGRAM CubeNBlk
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
          call Allgather_cube(x, blocksize, y, MPI_COMM_WORLD  )
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
      integer function log_base2( p) 
C  Just counts number of bits to right of most significant
C  bit.  So for p not a power of 2, it returns the floor
C  of log_2(p).
C  
      integer  p
      integer  return_val 
      integer  q
C
      return_val = 0
      q =  p
      do while(q .NE. 1) 
          q = q / 2
          return_val = return_val + 1
      end do
      log_base2 = return_val
      return
      end   
C
C
C ******************************************************************  
      subroutine Allgather_cube(x, blocksize, y, comm)
      real    x(0:blocksize -1)
      integer blocksize   
      real    y(0:1023)         
      integer comm
C
      include   'mpif.h'
      integer   i, d, p, my_rank
      integer   eor_bit
      integer   and_bits
      integer   stage, partner
      integer   hole_type
      integer   send_offset, recv_offset
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      integer   send_request
      integer   recv_request
      integer   commval, stride
C functions
      integer       log_base2
      integer       And
      integer       Xor
C
      call MPI_COMM_SIZE(comm,  p, ierr)
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
C  Copy x into correct location in y   
      do 100 i = 0 , blocksize  
          y(i + my_rank*blocksize) = x(i)
 100  continue
C
C  Set up   
      d = log_base2(p)
      eor_bit = 1
      and_bits = 1
      do 200 i = 1, d-1
           eor_bit = 2*eor_bit
 200  continue
      do 300 i = 1, d
           and_bits = 2*and_bits
 300  continue
      and_bits = and_bits - 1
C
      partner = XOr(my_rank ,eor_bit)
      send_offset = And(my_rank , and_bits)*blocksize
      recv_offset = And(partner ,and_bits)*blocksize
      call MPI_TYPE_CONTIGUOUS(blocksize, MPI_REAL, 
     +                    hole_type, ierr)
      call MPI_TYPE_COMMIT( hole_type, ierr)
C
      do 400  stage = 0 , d-1
          call MPI_ISEND(y(send_offset), 1, hole_type,
     +         partner, 0, comm,  send_request, ierr)
          call MPI_IRECV(y(recv_offset), 1, hole_type,
     +         partner, 0, comm,  recv_request, ierr)
C
          if (stage .LT. d-1)  then
              eor_bit = eor_bit/2
              and_bits = and_bits/2
              partner = XOr(my_rank ,eor_bit)
              send_offset = And(my_rank ,and_bits)*blocksize
              recv_offset = And(partner ,and_bits)*blocksize
              call MPI_TYPE_FREE( hole_type, ierr)
C
              commval = 1
              stride = 1
          do 500 i = 1,stage+1
              commval = commval * 2
 500      continue  
          do 600 i = 1, d-stage-1
              stride = stride * 2
 600      continue       

              call MPI_TYPE_VECTOR(commval, blocksize,
     +              stride*blocksize, MPI_REAL ,
     +              hole_type, ierr)
              call MPI_TYPE_COMMIT( hole_type, ierr)
          endif
C
          call MPI_WAIT( send_request,  status, ierr)
          call MPI_WAIT( recv_request,  status, ierr)
 400  continue
      return
      end
