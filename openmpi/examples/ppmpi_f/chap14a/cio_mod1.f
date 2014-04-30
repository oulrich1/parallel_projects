C  cio_mod1.f -- C  Collective functions for basic I/O
C      operations.  Modififed for use with parallel_bitonic.f
C
C  See Chap 8, pp. 142 & ff in PPMPI
C
C *******************************************************  
C  cio.c
C  Collective functions for basic I/O operations
C
C ******************************************************  
C  Attempt to identify a process in io_comm that can be
C      used for I/O.
C 
C  First see whether program defined io rank has been
C      cached with either communicator.  If this fails
C      try MPI defined io rank.
C 
C  Return values:
C      1.  0: rank of I/O process cached with io_comm.
C      2.  NO_IO_ATTR: couldn't find processor that could
C          carry out I/O.  MPI_PROC_NULL cached with
C          io_comm.
C 
C  Notes:
C      1.  This is a collective operation, since function
C          Copy_attr may use collective comm.
C      2.  Only possible values cached are a valid process
C          rank in comm2 or MPI_PROC_NULL.  (MPI_ANY_SOURCE
C          won't be cached
C
C*************************************************************
      integer function Cache_io_rank(orig_comm, io_comm)
      integer   orig_comm 
      integer   io_comm            
C
      include 'mpif.h'
      include 'ciof.h'
      integer retval
      integer my_rank
      integer ierr
C
      if (DEBUG) then
          call MPI_COMM_RANK(orig_comm,  my_rank, ierr)
          print *,'Process ' ,my_rank,' In Cache_io_rank'
      endif
C
C  Check whether IO_KEY is defined.  If not, define   
      if (IO_KEY .EQ. MPI_KEYVAL_INVALID)  then
          call MPI_KEYVAL_CREATE(MPI_DUP_FN,
     +          MPI_NULL_DELETE_FN,  IO_KEY,
     +          extra_arg, ierr)
      else
           retval = Copy_attr(io_comm, io_comm,IO_KEY)
           if (retval .NE. NO_IO_ATTR) then
C  Value cached   
              Cache_io_rank = retval
              return
            else 
               retval=Copy_attr(orig_comm,
     +                           io_comm, IO_KEY)
               if ( retval .NE. NO_IO_ATTR) then
C  Value cached   
                   Cache_io_rank = retval
                   return
               endif
            endif
      endif
C
      if (DEBUG) then
           print *, 'Process ',my_rank, 
     +               ' In Cache_io_rank, give up on IO_KEY'
      endif
C  Now see if we can find a value cached for MPI_IO 
      retval = Copy_attr(orig_comm, io_comm, MPI_IO)  
      if (retval .NE. NO_IO_ATTR) then
C  Value cached   
          if (DEBUG) then
               print *,'Process ',my_rank, 
     +                 ' > Copied attribute from orig to io'
                  
          endif
          Cache_io_rank =retval
          return
      else 
          retval = Copy_attr(io_comm, io_comm, MPI_IO)
          if (retval  .NE. NO_IO_ATTR) then
C  Value cached   
              if ( DEBUG) then
                 print *, 'Process ',my_rank, 
     +                  ' > Copied attribute from io to io'
              endif
              Cache_io_rank = retval
              return
          endif
      endif
C
      if ( DEBUG) then
            print *,'Process ',my_rank,
     +         ' > In Cache_io_rank, returning at end'
      endif
C  Couldn't find process that could carry out I/O   
C  Copy_attr has cached MPI_PROC_NULL 
      Cache_io_rank = NO_IO_ATTR              
      return 
      end   
C
C
C ******************************************************
C  Get attribute value associated with attribute key KEY
C      in comm1, and cache with comm2 IO_KEY
C 
C  KEY can be either IO_KEY or MPI_IO.
C 
C  Return values:
C      1.  0:  valid attribute successfully cached.
C      2.  NO_IO_ATTR:  Couldn't find process that could
C          carry out I/O.  MPI_PROC_NULL is cached with
C          comm2.
C 
      integer function Copy_attr(comm1, comm2, KEY)
      integer   comm1    
      integer   comm2    
      integer   KEY      
C
      include   'mpif.h'
      include   'ciof.h'
      integer   ierr 
      integer   io_rank
      integer   temp_rank
      integer   io_rank
      integer   equal_comm
      integer   flag 
      integer   my_rank
      integer   temp
C
      if (DEBUG) then
          call MPI_COMM_RANK(comm1,  my_rank, ierr)
          print *, 'Process ',my_rank,' > In Copy_attr'
      endif
C
      call MPI_ATTR_GET(comm1,KEY,io_rank,flag,ierr)
C
      if (DEBUG) then
          if (flag .EQ. 0)  then
              print *, 'Process ',my_rank,
     +                ' > Attr_get failed'
          else 
             if ( io_rank .EQ. MPI_ANY_SOURCE) then
                print *, 'Process ',my_rank,
     +                ' > attr = MPI_ANY_SOURCE'
                print *,'Process ',my_rank,
     +                           ' > MPI_ANY_SOURCE= ',
     +                           MPI_ANY_SOURCE
             else  
                print *,'Process ',my_rank,'> attr = ',
     +                    io_rank
             endif
         endif
      endif
C
      if (flag .EQ. 0)  then
          io_rank = MPI_PROC_NULL
          call MPI_ATTR_PUT(comm2, IO_KEY, 
     +       io_rank, ierr)
          Copy_attr = NO_IO_ATTR
          return
      else 
          if (io_rank .EQ. MPI_PROC_NULL)  then
              call MPI_ATTR_PUT(comm2, IO_KEY, 
     +             io_rank, ierr)
              Copy_attr = NO_IO_ATTR
              return
          else 
              if (io_rank .EQ. MPI_ANY_SOURCE) then
C  Any process can carry out I/O.  Use   
C  process 0                             
                 if( DEBUG) then
                     print *, 'Process ',my_rank,
     +              ' > Return from Copy on MPI_ANY_SOURCE'
                 endif
                 io_rank = 0
                 call MPI_ATTR_PUT(comm2, IO_KEY, 
     +                       io_rank, ierr)
                 if ( DEBUG) then
                     call MPI_ATTR_GET(comm2, IO_KEY,  
     +                 temp,  flag, ierr)
                     if (flag .EQ. 0) then
                           print *, 'Process ',my_rank,
     +                     ' > In Copy, no value cached!'
                     else
                           print *, 'Process ',my_rank,
     +                     ' > In Copy, cached io_rank = ',
     +                    temp
                     endif
                 endif
                 Copy_attr = 0
                 return
              endif
          endif
      endif
C
C  Value in *io_rank is a valid process    
C  rank in comm1.  Action depends on whether   
C  comm1 == comm2.                             
      call MPI_COMM_COMPARE(comm1, comm2,  
     +             equal_comm, ierr)
C
      if (equal_comm .EQ. MPI_IDENT)  then
C  comm1 == comm2.  Valid value already   
C  cached.  Do nothing.                   
          Copy_attr = 0
      else  
C  Check whether rank returned is valid   
C  process rank in comm2                  
          call Get_corresp_rank(comm1, io_rank,
     +         comm2,  temp_rank)
C
C  Different ranks may have been returned   
C  on different processes.  Get min         
          if (temp_rank .EQ. MPI_UNDEFINED) then
              temp_rank = HUGE
          endif
          print*, 'Process ',my_rank,'> temp_rank = ', temp_rank
          call MPI_ALLREDUCE( temp_rank,  io_rank, 1,
     +         MPI_INTEGER, MPI_MIN, comm2, ierr)
C
          if (io_rank .LT. HUGE)  then
              io_rank = io_rank
              call MPI_ATTR_PUT(comm2, IO_KEY, 
     +                 io_rank, ierr)
              Copy_attr = 0
          else  
C  No process got a valid rank in comm2   
C  from Get_corresp_rank                  
              io_rank = MPI_PROC_NULL
              call MPI_ATTR_PUT(comm2, IO_KEY, 
     +                 io_rank,ierr)
              Copy_attr = NO_IO_ATTR               
          endif
      endif
      end   
C
C ******************************************************  
C  Determines whether the process with rank rank1 in 
C      comm1 is a valid rank in comm2.
C  If it is, it returns the rank in *rank2.  If it 
C      isn't it returns MPI_UNDEFINED.
C 
      subroutine Get_corresp_rank(comm1, rank1,
     +                        comm2, rank2)
      integer   comm1        
      integer   rank1        
      integer   comm2        
      integer   rank2    
C
      include  'mpif.h'
      integer   ierr 
      integer   group1, group2
C
      call MPI_COMM_GROUP(comm1,  group1, ierr)
      call MPI_COMM_GROUP(comm2,  group2, ierr)
C
      call MPI_GROUP_TRANSLATE_RANKS(group1, 1,
     +          rank1, group2, rank2, ierr)
C
      return
      end   
C
C ******************************************************  
C  Check whether IO_KEY is valid.  If it is, attempt to
C      access it.  If it isn't attempt to define it from
C      MPI_COMM_WORLD.
C 
C  Return values:
C      1.  0: Valid rank returned.
C      2.  NO_IO_ATTR:  Unable to find rank.
C 
      integer function Get_io_rank(io_comm, io_rank)
      integer io_comm       
      integer io_rank   
C
      include  'mpif.h'
      include  'ciof.h'
      integer  ierr 
      integer  temp
      integer  flag
      integer  retval
C
      if (IO_KEY .EQ. MPI_KEYVAL_INVALID)  then
          call MPI_KEYVAL_CREATE(MPI_DUP_FN,
     +         MPI_NULL_DELETE_FN,  IO_KEY,
     +         extra_arg, ierr)
      else  
          call MPI_ATTR_GET(io_comm, IO_KEY,temp,  
     +           flag, ierr)
          if ((flag .NE. 0) .AND.
     +        (temp .NE. MPI_PROC_NULL))  then
              io_rank = temp
              Get_io_rank = 0
              return
          endif
      endif
C
      retval = Copy_attr(MPI_COMM_WORLD, io_comm,
     +                       MPI_IO)
      if (retval .EQ. NO_IO_ATTR) then
          Get_io_rank = NO_IO_ATTR
          return
      else  
          call MPI_ATTR_GET(io_comm, IO_KEY, temp,
     +         flag, ierr)
          io_rank = temp
          Get_io_rank = 0
          return
      endif
C
      end   
C
C ******************************************************  
C  TEMPLATE - see comments for needed changes!
C  Prompt for input, read one line and broadcast.  
C 
C  Return values:
C      1. 0:  input read
C      2. NO_IO_ATTR:  no rank cached with IO_KEY.
C 
C  Notes:
C      1. Prompt is significant only on IO_process 
C 
C      User adds paramters between 2 & N as needed
C
C Template for reading in data.      
      integer function Cread(io_comm, N, prompt, 
     +                    param1 )     
      integer      io_comm
      integer      N
      character *40 prompt
C DEFAULT type for parameters. User updates as needed!
      integer      param1
      include      'mpif.h'
      include      'ciof.h'
      integer      my_io_rank
      integer      root                                  
      integer      derived_mpi_t
      integer      ierr
      integer      retval
C
      call MPI_COMM_RANK(io_comm, my_io_rank, ierr)
C
      retval = Get_io_rank(io_comm, root)
      if (retval .eq. NO_IO_ATTR) then
          Cread =  NO_IO_ATTR
          return
      endif
      if (my_io_rank .EQ.  root) then
          print *, prompt
          read *,  param1 
      endif
      call MPI_BCAST(param1, 1, MPI_INTEGER, root,
     +                  io_comm, ierr)
      Cread = 0   
      return
      end
C 
Cread
C***********************************************************************
C  TEMPLATE ONLY!! 
C  Here user is also required to add paramaters as needed and
C  update the types (default is integer)
C  
      subroutine Build_derived_type(N, param1,  param2, paramN,
     +     derived_mpi_t)
      integer   HUGE
      parameter (HUGE = 20)
      integer   N
C default types - update as needed
      integer   param1
      real      param2
C add additional parameters here. .between 2 & N
      integer   paramN
C
      include   'mpif.h'
      integer   derived_mpi_t
      integer   block_lengths(HUGE)
      integer   displacements(HUGE)
      integer   types(HUGE)
      integer   ierr      
      integer   i
      integer   start_address
      integer   curr_address
C
      do 100 i = 1, N
          block_lengths(i) = 1
  100 continue
C
C Change MPI_types to match the type for each parameter
      types(1) = MPI_INTEGER
      types(2) = MPI_REAL
      types(3) = MPI_INTEGER
C      
      call MPI_ADDRESS(param1, start_address, ierr)
      displacements(1) = 0
      call MPI_ADDRESS(param2, curr_address, ierr)
      displacements(2) = curr_address - start_address
C Repeat logic for any additional paramters
      call MPI_ADDRESS(paramN, curr_address, ierr)
      displacements(3)=curr_address - start_address
      call MPI_TYPE_STRUCT(N, block_lengths,
     +      displacements,
     +      types, derived_mpi_t, ierr)
      call MPI_TYPE_COMMIT(derived_mpi_t, ierr)
C
      return
      end
C Build_derived_type
C
C ******************************************************  
C  TEMPLATE ONLY!! USer updates for their environment!
C  Prints data from all processes.  Format of data must
C      be the same on each process.
C 
C  Return values:
C      1.  0:  data printed
C      2.  NO_IO_ATTR:  no rank cached with IO_KEY
C 
C  Notes:
C      1.  Title is significant only on root.
C 
C user would add paramaters here...
      integer function  Cprint(io_comm, num_params, title, 
     +                    param1 )
      integer  io_comm   
      integer        num_params
      character *40  title     
C types & num_params updated by user as appropriate
      integer        param1(128)
C
      include 'ciof.h'
      include 'mpif.h'
      integer         q
      integer         my_io_rank
      integer         io_p
      integer         root
      integer         status(MPI_STATUS_SIZE)
      integer         ierr
      integer         retval
      integer         derived_mpi_type
      integer         i
C
      call MPI_COMM_RANK(io_comm,  my_io_rank, ierr)
C
      retval = Get_io_rank(io_comm,  root)
      if ( retval .EQ. NO_IO_ATTR) then
          Cprint = NO_IO_ATTR
          return
      endif
      call MPI_COMM_SIZE(io_comm,  io_p, ierr)
C
C  Send output data to io_process   
      if (my_io_rank .NE. root)  then
C  Copy the output data into io_buf
          call MPI_SEND(param1, num_params ,
     +            MPI_INTEGER,
     +            root, 0, io_comm, ierr)              
    
      else    
          do 100 q =0, root-1
             call MPI_RECV(param1, num_params ,
     +                MPI_INTEGER, 
     +                q, 0, io_comm,  status, ierr)
             print *, ' ',title  
             print *,'Process ',q ,' > '
             print *, (param1(i),i=1,num_params)
     +                 
 100      continue
     
C
C  Copy the output data into io_buf
          print *, 'Process',root, ': ',title 
          print *,(param1(i),i=1,num_params)
C        
           do 300 q = root+1, io_p -1   
               call MPI_RECV(param1, num_params,
     +                 MPI_INTEGER,
     +                 q, 0, io_comm,  status, ierr)
               print *,'Process ',q, ' > ' 
               print *,(param1(i),i=1,num_params)
 300       continue
 
           print *, ' '
      endif
C  
      Cprint = 0
      return 
      end  
C
C
C ******************************************************  
C  Gathers error codes from all processes to all processes.
C      If any error is negative, all processes abort.   
C 
C  Return values:
C      1. 0:  no error detected.
C      2. NO_IO_ATTR:  No valid rank cached with IO_KEY.
C 
C  Notes:
C      1. "routine_name" only has significance on io_process. 
C
      integer function Cerror_test(io_comm, 
     +                       routine_name, error)
      integer     io_comm        
      character   *30  routine_name   
      integer     error          
C
      include 'mpif.h'
      include 'ciof.h'
      integer q
      integer io_p
      integer error_count
      integer io_process
      integer my_io_rank
      integer error_buff(0:99)
      integer ierr
C
      error_count = 0
      if (Get_io_rank(io_comm,  io_process) .EQ. 
     +               NO_IO_ATTR)   then
         Cerror_test = NO_IO_ATTR
         return
      endif
      call MPI_COMM_SIZE(io_comm,  io_p, ierr)
      call MPI_COMM_RANK(io_comm,  my_io_rank, ierr)
C
      call MPI_ALLGATHER( error, 1, MPI_INTEGER, error_buf, 
     +     1, MPI_INTEGER, io_comm, ierr)
C
      do 100 q = 0, io_p-1 
          if (error_buff(q) .LT. 0)  then
              error_count = error_count +1
              if (my_io_rank .EQ. io_process)  then
                  print *, 'Error in ', routine_name,
     +             ' on process ', q
              endif
          endif
 100  continue
C
      if (error_count .GT. 0) then
          call MPI_ABORT(MPI_COMM_WORLD, -1, ierr )
      endif
C
      Cerror_test = 0
      return
      end

