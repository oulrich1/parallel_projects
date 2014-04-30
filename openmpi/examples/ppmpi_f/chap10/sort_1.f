C  sort_1.f -- level 1 version of sort program
C 
C  Input: none
C  Output:  messages indicating flow of control through program
C 
C  See Chap 10, pp. 226 & ff in PPMPI.
C 
      PROGRAM Sort1
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
      INCLUDE 'sort_1f.h'
      integer       local_keys(0:1000)
      integer       listsize
      integer       error
      integer       ierr
      integer       retval
      integer       p
      integer       my_rank
      integer       io_comm
C
       IO_KEY = MPI_KEYVAL_INVALID    
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      retval =  Cache_io_rank(MPI_COMM_WORLD, io_comm)
C
      listsize = Get_list_size( my_rank)
C
C  Return negative if Allocate failed   
      error = Allocate_list(listsize, local_keys, my_rank)
C
      call Get_local_keys( local_keys,my_rank)
      call Print_list( local_keys ,my_rank)
      call Redistribute_keys( local_keys,my_rank)
      call Local_sort( local_keys,my_rank )
      call Print_list( local_keys,my_rank )
C
      call MPI_FINALIZE(ierr)
      end
C******************************************************************
      integer function Get_list_size( my_rank)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Get_list_size'
      Get_list_size = 0
      return 
      end
C
C  Return value negative indicates failure   
      integer function Allocate_list(list_size, local_keys,
     +           my_rank )
      integer list_size
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Allocate_key_list'
      Allocate_list = 0
      return 
      end  
C
      subroutine Get_local_keys(  local_keys, my_rank )
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Get_local_keys' 
      return
      end  
C
      subroutine Redistribute_keys( local_keys, my_rank)
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Redistribute_keys' 
      return
      end  
C
      subroutine Local_sort(  local_keys, my_rank)  
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Local_sort'
      return
      end  
C
      subroutine Print_list(  local_keys , my_rank)  
      integer local_keys(0:1000)
      integer       my_rank
      print *, 'Proc ' ,my_rank, '> In Print_list' 
      return
      end
