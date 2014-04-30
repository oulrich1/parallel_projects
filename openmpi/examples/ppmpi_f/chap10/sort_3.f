C  sort_3.f -- level 3 version of sort program
C  2.  Add fake Get_list_size 
C  2.  Return scalars in Allocate_list
C  2.  Add fake Get_local_keys
C  3.  Definition of Redistribute keys.
C  3.  Finish Allocate_list.
C  3.  Add Insert.
C  3.  Add Local_sort.
C  3.  Add Print_list.
C 
C  Input: none
C  Output: messages indicating flow of control and contents of
C      local lists.
C 
C  Note:  Doesn't call Redistribute_keys.
C 
C  See Chap 10, pp. 226 & ff, esp. pp. 232 & ff., in PPMPI.
C
C *******************************************************************  
      PROGRAM Sort3
      INCLUDE 'mpif.h'
      INCLUDE 'sort_3f.h'
      INCLUDE 'ciof.h'
      integer  local_keys(0:1000)
      data local_keys /1001 * 0/
      integer       listsize
      integer       error
      integer       ierr
      integer       retval
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
      listsize = Get_list_size(p)
C
C  Return negative if Allocate failed   
      error = Allocate_list(listsize,  local_keys, p)
      if (error .LT. 0)  then
          print *, 'Process ',my_rank,' > Can''t allocate list ' 
          print *, 'Process ',my_rank,' > Quitting '
          call MPI_ABORT(MPI_COMM_WORLD, -1, ierr )
      endif
C
      call Get_local_keys( local_keys, my_rank)
      call Print_list(io_comm,  local_keys)
C 
      call Redistribute_keys( local_keys, p)
      call Local_sort( local_keys(list_size),local_keys)
      call Print_list(io_comm,  local_keys)
C
      call MPI_FINALIZE(ierr)
      end
C
C
C******************************************************************
      integer function Get_list_size(p )
      integer       p
      integer       my_rank
      integer       io_comm
      common        my_rank, io_comm
      print *, 'Proc ' ,my_rank, '> In Get_list_size'
      Get_list_size = 5*p
      return 
      end
C******************************************************************
C  Return value negative indicates failure   
      integer function Allocate_list(listsize, local_keys, p )
      integer listsize
      integer       p
      integer local_keys(0:1000)
      include       'sort_3f.h'
C
      local_keys(allocated_size) = listsize/p
      local_keys(list_size) = listsize/p
      Allocate_list = 0
      return 
      end  
C
C******************************************************************
      subroutine Get_local_keys(  local_keys, my_rank )
      integer local_keys(0:1000)
      include 'sort_3f.h'
      integer       my_rank
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
C *******************************************************************  
      subroutine Redistribute_keys(local_keys, p)
      integer local_keys(0:1000)
      integer p
      include 'sort_3f.h'
      integer new_list_size, i, error 
      data error /0/
      include 'mpif.h'
      integer ierr
C 63 below assumes max processors = 64.  change if needed
      integer send_counts(0:63)
      integer send_displacements(0:63)
      integer recv_counts(0:63)
      integer recv_displacements(0:63)
      integer new_keys(0:1000)
      data new_keys /1001*0/
C
      call Local_sort(local_keys(list_size), local_keys(keys))
      call Find_alltoall_send_params(local_keys,
     +     send_counts, send_displacements)
C
C  Distribute the counts   
      call MPI_ALLTOALL(send_counts, 1, MPI_INTEGER,
     +     recv_counts,
     +     1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
C
C  Allocate space for new list   
      new_list_size = recv_counts(0)
      do 100 i = 1, p-1
          new_list_size = new_list_size + recv_counts(i)
 100  continue
C
       call Find_recv_displacements(recv_counts, 
     +                             recv_displacements)
C
C  Exchange the keys   
       call MPI_ALLTOALLV(local_keys(keys), send_counts,
     +     send_displacements, MPI_INTEGER, new_keys,
     +     recv_counts, recv_displacements, MPI_INTEGER,
     +     MPI_COMM_WORLD, ierr )
C
C  Replace old list with new list   
       local_keys(allocated_size) = new_list_size
       local_keys(list_size) = new_list_size
       do 200 i = 0, new_list_size
          local_keys(keys+i) = new_keys(i)
 200   continue
C
       return
       end
C
C *******************************************************************  
      subroutine Find_alltoall_send_params(local_keys,
     +     send_counts, send_displacements)
      integer  local_keys(0:1000)
C 63 below assumes max processors = 64.  change if needed
      integer send_counts(0:63)
      integer send_displacements(0:63)
C do nothing...set to zero all values
      integer i
      do 100 i = 0, 63
         send_counts(i) = 0
         send_displacements(i)= 0
 100  continue
      return
      end  
C
C *******************************************************************  
      subroutine Find_recv_displacements(recv_counts,
     +                recv_displacements ) 
C 63 below assumes max processors = 64.  change if needed
      integer recv_counts(0:63)
      integer recv_displacements(0:63)
C do nothing... 
      return
      end
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
      if ( p.LT. q) then
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
C *******************************************************************  
      subroutine Print_list(io_comm,  local_keys) 
      integer io_comm
      integer local_keys(0:1000)
      include 'ciof.h'
      include 'sort_3f.h'
      integer  i
      character *40 title
      integer retval
      title = 'Contents of the list                      '
      retval = Cprint(io_comm, local_keys(list_size),
     +             title, local_keys(keys))
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
