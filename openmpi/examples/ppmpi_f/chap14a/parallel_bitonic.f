C  parallel_bitonic.f -- parallel bitonic sort of randomly generated list
C      of integers
C 
C  Input:
C      n: the global length of the list -- must be a power of 2.
C 
C  Output:
C      The sorted list.
C 
C  Notes:
C      1.  Assumes the number of processes p = 2^d and p divides n.
C      2.  The lists are statically allocated -- size specified in MAX.
C      3.  Keys are in the range 0 -- KEY_MAX-1.
C      4.  Implementation can be made much more efficient by using
C          pointers and avoiding re-copying lists in merges.
C 
C  See Chap 14, pp. 320 & ff. in PPMPI.
C
C ******************************************************************
      PROGRAM ParBitonic
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
C
       integer MAX  
       integer LOW
       integer HIGH
       integer KEY_T
       integer KEY_MAX
       parameter (MAX = 16384, low = 0, high = 1, KEY_MAX = 32768)
C
      integer       list_size           
      integer       n                   
      integer       local_list(0:16383)
      integer       proc_set_size
      integer       my_rank
      integer       p
      integer       and_bit
      integer       io_comm
      integer       ierr
      integer       retval
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval = Cache_io_rank(MPI_COMM_WORLD, io_comm  )
C
      retval = 
     +   Cread(io_comm,1, 'Enter the global list size         ', n)
      list_size = n/p
C
      call Generate_local_list(list_size, local_list)
C 
      call Print_list('Before local sort ', list_size, 
     +                local_list, io_comm)
C        
      call Local_sort(list_size, local_list)
C 
      call Print_list('After local sort  ', list_size, 
     +                local_list, io_comm)
C 
C
C  and_bit is a bitmask that, when "anded" with    
C  my_rank, tells us whether we're working on an   
C  increasing or decreasing list                   
      proc_set_size = 2
      and_bit = 2
      do while ( proc_set_size .LE. p)
          if (AND(my_rank ,and_bit) .EQ. 0) then
              call Par_bitonic_sort_incr(list_size,
     +                   local_list, proc_set_size, 
     +                   MPI_COMM_WORLD  )
          else
              call Par_bitonic_sort_decr(list_size,
     +                   local_list, proc_set_size, 
     +                   MPI_COMM_WORLD )
          endif
          proc_set_size = proc_set_size * 2
          and_bit = and_bit * 2
      end do
C
      call Print_list('After sort        ', list_size, local_list, 
     +                     io_comm)
C
      call MPI_FINALIZE(ierr)
      end
C
C
C *******************************************************************  
      subroutine Generate_local_list(list_size, local_list)
      integer    list_size      
      integer    local_list(0:list_size-1)   
      include    'mpif.h'
      integer i
      integer my_rank
      integer ierr
      integer iseed
      integer retval
      integer KEY_T
      integer KEY_MAX, MAX
      parameter (MAX = 16384,  KEY_MAX = 32768)
C function
      integer Random
C
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      iseed = 5 + my_rank
C
      do 100 i = 0 , list_size -1
          retval = Random(  iseed)
          local_list(i) = MOD (retval, KEY_MAX)
 100  continue 
C
      return
      end
C
C *******************************************************************  
      subroutine Print_list( title, list_size, local_list, io_comm)
      character *18 title         
      integer       list_size     
      integer       local_list(0:list_size-1)  
      integer       io_comm 
C
      include         'mpif.h'
      include         'ciof.h'     
      integer         i, q
      integer         p
      integer         my_rank
      integer         root
      integer         status(MPI_STATUS_SIZE)
      integer         ierr
      integer         retval
C
       integer temp_list(0:16383)
 
C
      call MPI_COMM_SIZE(io_comm,  p, ierr)
      call MPI_COMM_RANK(io_comm,  my_rank, ierr)
      retval = Get_io_rank(io_comm,  root)
C
      if (my_rank .EQ. root)  then
          print *, title 
          do 100 q = 0 ,root-1
              call MPI_RECV(temp_list, list_size, MPI_INTEGER,
     +            q, 0, io_comm,  status, ierr)
              print *,'Process ',q , '>'
              print *, (temp_list(i), i = 0,list_size-1)
 100      continue
          print *, 'Process ',root,' > '
          print *,( local_list(i), i = 0,list_size-1)
          do 200 q = root+1, p-1
              call MPI_RECV(temp_list, list_size, MPI_INTEGER,
     +             q, 0, io_comm,  status, ierr)
              print *, 'Process ',q ,' > '
              print *,( temp_list(i), i = 0,list_size-1)
 200      continue
      else  
          call MPI_SEND(local_list, list_size, MPI_INTEGER,
     +         root, 0, io_comm, ierr)
      endif
C
      return
      end   
C
C
C *******************************************************************  
      subroutine Local_sort(list_size, local_keys)
C   inefficient, but easy to code, local bubble sort
      integer    list_size      
      integer    local_keys( list_size  ) 
      integer    i, pass
      integer    temp
      logical    sorted  
C
      sorted = .FALSE.
      pass = 1
      do while (.NOT. sorted)
          sorted = .TRUE.
          do 100 i = 1, list_size - pass
              if (local_keys(i)
     +                   .GT. local_keys(i+1)) then
C swap
                 temp = local_keys(i)
                 local_keys(i) = local_keys(i + 1)
                 local_keys(i + 1) = temp
                 sorted = .FALSE.
              endif
 100      continue
          pass = pass + 1
      end do
C
      return
      end        
C *******************************************************************  
      integer function Key_compare(  p,   q) 
      integer p, q
C
      if ( p .LT. q) then
          Key_compare = -1
      else 
          if (p .EQ. q) then
            Key_compare = 0
          else  
            Key_compare = 1
          endif
      endif
C
      return
      end   
C ******************************************************************  
      integer function log_base2(x) 
      integer x
      integer xinit
      integer count 
C
      count = 0
      xinit = x
      do while (xinit .GT. 1) 
          xinit = xinit/2
          count = count + 1
      end do
C
      log_base2 = count
      return
      end
C   
C ******************************************************************  
      subroutine Par_bitonic_sort_incr(list_size, local_list,
     +                  proc_set_size, comm)
      integer       list_size       
      integer       local_list(0:list_size-1)     
      integer       proc_set_size   
      integer       comm            
C
      include       'mpif.h'
      integer LOW
      integer HIGH
      parameter (  low = 0, high = 1)
      integer       eor_bit
      integer       proc_set_dim
      integer       stage
      integer       partner
      integer       my_rank
      integer       ierr
      integer       i
C functions
      integer       log_base2
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      proc_set_dim = log_base2(proc_set_size)
      eor_bit = 1
      do 100 i = 1, proc_set_dim -1
         eor_bit = eor_bit * 2
 100  continue
      do 200 stage = 0 ,proc_set_dim  -1
          partner = XOr(my_rank, eor_bit)
          if (my_rank .LT. partner) then
              call Merge_split(list_size, local_list, LOW,
     +             partner, comm)
          else
              call Merge_split(list_size, local_list, HIGH,
     +             partner, comm)
          endif
          eor_bit = eor_bit/2
 200  continue
C
      return
      end  
C
C
C
C ******************************************************************  
      subroutine Par_bitonic_sort_decr(list_size, local_list,
     +                  proc_set_size, comm)
      integer       list_size       
      integer       local_list(0:list_size-1)      
      integer       proc_set_size   
      integer       comm            
C
      include       'mpif.h'
      integer LOW
      integer HIGH
      parameter (  low = 0, high = 1)
      integer       eor_bit
      integer       proc_set_dim
      integer       stage
      integer       partner
      integer       my_rank
      integer       ierr
      integer       i
C functions
      integer       log_base2
C
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      proc_set_dim = log_base2(proc_set_size)
      eor_bit = 1
      do 100 i = 1, proc_set_dim -1
         eor_bit = eor_bit * 2
 100  continue
      do 200 stage = 0 ,proc_set_dim  -1
          partner = XOr(my_rank, eor_bit)
          if (my_rank .GT. partner) then
              call Merge_split(list_size, local_list, LOW,
     +             partner, comm)
          else
              call Merge_split(list_size, local_list, HIGH,
     +             partner, comm)
          endif
          eor_bit = eor_bit/2
 200  continue
C
      return
      end  
C
C
C ******************************************************************  
      subroutine Merge_split(list_size, local_list, which_keys,
     +    partner, comm)
      integer       list_size      
      integer       local_list(0:list_size-1)   
      integer       which_keys
      include       'mpif.h'    
      integer       partner        
      integer       comm   
      integer LOW
      integer HIGH
      parameter ( low = 0, high = 1 )
C
       integer temp_list(0:16383)
C
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
C
C  key_mpi_t is an MPI (derived) type   
      call MPI_SENDRECV(local_list, list_size, MPI_INTEGER,
     +          partner, 0, temp_list, list_size,
     +          MPI_INTEGER, partner, 0, comm,  status, ierr)
      if (which_keys .EQ. HIGH) then
          call Merge_list_high(list_size, local_list, temp_list)
      else
          call Merge_list_low(list_size, local_list, temp_list)
      endif
      return
      end  
C
C ******************************************************************  
C  Merges the contents of the two lists.   
C  Returns the smaller keys in list1       
      subroutine Merge_list_low(list_size, list1, list2)
      integer    list_size   
      integer    list1(0:list_size-1)
      integer    list2(0:list_size-1)
      integer  i
      integer  index1  
      integer  index2  
      data index1, index2 / 0, 0/
C  in Merge_split
       integer scratch_list(0:16383)
       common scratch_list
C
      do 100 i = 0  , list_size-1
          if (list1(index1) .LE. list2(index2))  then
              scratch_list(i) = list1(index1)
              index1 = index1 + 1
          else  
              scratch_list(i) = list2(index2)
              index2 = index2 + 1
          endif
 100  continue
      do 200 i =0, list_size -1
          list1(i) = scratch_list(i)
 200  continue
C
      return 
      end   
C
C
C ******************************************************************  
C  Returns the larger keys in list 1.      
      subroutine Merge_list_high(list_size, list1, list2)
      integer    list_size
      integer    list1(0:list_size-1)
      integer    list2(0:list_size-1)
      integer  i
      integer  index1  
      integer  index2 
C  in Merge_split
       integer scratch_list(0:16383)
       common scratch_list
C
      index1 = list_size -1
      index2 = list_size -1
      do 100 i = list_size - 1, 0 , -1
          if (list1(index1) .GE.  list2(index2))  then
              scratch_list(i) = list1(index1)
              index1 = index1 - 1
          else  
              scratch_list(i) = list2(index2)
              index2 = index2 - 1
          endif
 100  continue
      do 200 i = 0, list_size-1
          list1(i) = scratch_list(i)
 200  continue
C
      return
      end   
C
C**********************************************************
C
      integer function Random(seed)
C
      integer seed
      integer c1, c2
      parameter (c1 = 19423, c2 = 3011)
      integer oldseed
      save oldseed
      data oldseed /0/
C
      if (oldseed .EQ. 0) then
          oldseed = seed
      endif
C
      oldseed = mod(c1*oldseed, c2)
C
      Random = oldseed 
C
      return
      end
