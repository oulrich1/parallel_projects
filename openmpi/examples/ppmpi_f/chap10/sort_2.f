C  sort_2.f -- level 2 version of sort program
C  2. Add fake Get_list_size 
C  2. Return scalars in Allocate_list 
C  2. Add fake Get_local_keys
C 
C  Input: none
C  Output:  messages indicating flow of control
C 
C  See Chap 10, pp. 226 & ff, esp. pp. 231 & ff, in PPMPI
C
C *******************************************************************  
      PROGRAM Sort2
      INCLUDE 'mpif.h'
      INCLUDE 'sort_2f.h'
      INCLUDE 'ciof.h'
      integer  local_keys(0:1000)
      integer           listsize
      integer           error
      integer           ierr
      integer           retval
C
      integer       p
      integer       my_rank
      integer       io_comm
C
C
      IO_KEY = MPI_KEYVAL_INVALID    
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval =  Cache_io_rank(MPI_COMM_WORLD, io_comm)
C
      listsize = Get_list_size(my_rank , p)
C
C  Return negative if Allocate failed   
      error = Allocate_list(listsize,  local_keys,p )
C
      call Get_local_keys( local_keys )
      call Print_list( local_keys, my_rank )
      call Redistribute_keys( local_keys, my_rank )
      call Local_sort( local_keys, my_rank )
      call Print_list( local_keys, my_rank )
C
      call MPI_FINALIZE(ierr)
      end
C
C******************************************************************
      integer function Get_list_size(my_rank, p )
      integer       my_rank, p
      print *, 'Proc ' ,my_rank, '> In Get_list_size'
      Get_list_size = 5*p
      return 
      end
C******************************************************************
C  Return value negative indicates failure   
      integer function Allocate_list(listsize, local_keys,p )
      integer listsize
      integer local_keys(0:1000)
      include       'sort_2f.h'
      integer       p
C
      local_keys(allocated_size) = listsize/p
      local_keys(list_size) = listsize/p
      Allocate_list = 0
      return 
      end  
C
C******************************************************************
      subroutine Get_local_keys(  local_keys )
      integer local_keys(0:1000)
      include       'sort_2f.h'
      integer       p
      integer       my_rank
      integer       io_comm
      common        p, my_rank, io_comm
      integer       iseed, randval, i
C function
      integer       Random
C  Seed the generator   
      iseed = 5 + my_rank
      do 100 i = 0 , local_keys(list_size)-1    
          randval = Random(iseed)
          call Insert_key(MOD(randval,KEY_MOD),i,local_keys)
 100  continue 
C
      return
      end  
C
C******************************************************************
      subroutine Insert_key(key,  i, local_keys)
      integer key
      integer i
      integer local_keys(0:1000) 
C do nothing...
      return
      end
C******************************************************************
      subroutine Redistribute_keys( local_keys, my_rank)
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Redistribute_keys' 
      return
      end  
C
C******************************************************************
      subroutine Local_sort(  local_keys, my_rank)  
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Local_sort'
      return
      end  
C
C******************************************************************
      subroutine Print_list(  local_keys, my_rank )  
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Print_list' 
      return
      end 
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
