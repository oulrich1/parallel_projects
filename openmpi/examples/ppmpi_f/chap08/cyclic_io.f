C  cyclic_io.f -- Functions for I/O of arrays using a cyclic 
C      distribution.
C 
C  See Chap 8, pp. 158 & ff in PPMPI
C 
C ******************************************************  
C  Initialize all members except entries
C 
      subroutine Initialize_params(comm, n, array)
      integer      comm    
      integer      n       
      integer      array(0: 8 + 1024 + 1024-1)   
C
      INCLUDE   'mpif.h'
      INCLUDE   'cyclic_iof.h'
      integer   p
      integer   my_rank
      integer   q
      integer   quotient
      integer   remainder
      integer   ierr
C
      call MPI_COMM_SIZE(comm,  p, ierr)
      call MPI_COMM_RANK(comm,  my_rank, ierr)
C
      array(arr_comm) = comm
      array(arr_p) = p
      array(arr_myrank) = my_rank
C
      array(global_order) = n
C
      quotient = n/p
      remainder = MOD(n , p)
C
      if (remainder .EQ. 0) then
          array(padded_size) = n
      else
          array(padded_size) = p*(quotient+1)
      endif
C 
      if (my_rank .LT. remainder)  then
          array(local_size) = quotient+1
      else  
          array(local_size) = quotient
      endif
C
      array(arr_stride) = p
C
      call Build_cyclic_type(array(type), array(arr_stride),
     +                    array(padded_size), p)
C
      return
      end    
C
C
C ******************************************************  
      subroutine Build_cyclic_type(cyclic_mpi_t, stride,
     +                             array_size, p)
      integer           cyclic_mpi_t   
      integer           stride         
      integer           array_size     
      integer           p
      INCLUDE           'mpif.h'
      integer           ierr           
C
C change for your system size!
      integer           sizeofint
      parameter         (sizeofint = 16)
      integer           vector_mpi_t
      integer           blocksizes(0:1)
      integer           displacements(0:1)
      integer           type_list(0:1)
C
      call MPI_TYPE_VECTOR(array_size/p, 1, stride, 
     +                    MPI_INTEGER, vector_mpi_t, ierr)
C
      blocksizes(0) = 1
      blocksizes(1) = 1
      displacements(0) = 0
      displacements(1) = sizeofint
      type_list(0) = vector_mpi_t
      type_list(1) = MPI_UB
C
      call MPI_TYPE_STRUCT(2, blocksizes, displacements, 
     +                  type_list, cyclic_mpi_t, ierr)
      call MPI_TYPE_COMMIT(cyclic_mpi_t, ierr)
      return
      end 
C
C ******************************************************  
      subroutine Print_params(array)
      integer  array(0: 8 + 1024 + 1024-1)
      INCLUDE  'cyclic_iof.h'
      INCLUDE  'ciof.h'
      character *15 title
      integer   retval
C
      title = 'p =           '
      retval = Cprint(array(arr_comm), 1, title, array(arr_p))
      title = 'my_rank =     '
      retval = Cprint(array(arr_comm), 1, title, array(arr_myrank))
      title = 'order =       '
      retval = Cprint(array(arr_comm), 1, title, array(global_order))
      title = 'padded_size = '
      retval = Cprint(array(arr_comm), 1, title, array(padded_size))
      title = 'my_size =     '
      retval = Cprint(array(arr_comm), 1, title, array(local_size))
      title = 'stride =      '
      retval = Cprint(array(arr_comm), 1, title, array(arr_stride))
C
      return
      end  
C
C ******************************************************  
C  Assumes that each process is using local_entries 
C     member to store current contents of array.  If
C     this is not the case, appropriate range of
C     values from entries must be copied into 
C     local_entries before call to MPI_Gather.
C 
      subroutine Print_entries(title, array)
      character *(*)          title   
      integer  array(0: 8 + 1024 + 1024 -1)   
      INCLUDE 'cyclic_iof.h'
      INCLUDE 'ciof.h'
      INCLUDE 'mpif.h'
C
      integer root
      integer q
      integer quotient
      integer remainder
      integer i, j, k
      integer send_size
      integer incr
      integer   prtarr(0:63)
      character *12 prtlines(0:63)
      integer retval
      integer ierr
C
      retval = Get_io_rank(array(arr_comm),  root)
C
      send_size = array(padded_size)/array(arr_p)
      call MPI_GATHER(array(local_entry), send_size, MPI_INTEGER,
     +        array(entry), 1, array(type), root,
     +        array(arr_comm), ierr)
C
      if (array(arr_myrank) .EQ. root)  then
          print *, title
          print *, '   Processes '
          do 100 q = 0 , array(arr_p)-1   
              prtarr(q) = q
              prtlines(q) = '    --------'
 100      continue
          print *, (prtarr(q), q=0,array(arr_p) -1)
          print *, (prtlines(q), q=0,array(arr_p)-1)
C
          quotient = array(global_order)/array(arr_p)
          remainder = MOD( array(global_order) , array(arr_p) )
C
          do 400 i = 0 , quotient-1   
              print *, (array(entry + (i*array(arr_p)) +j ),
     +                    j = 0, array(arr_p)-1 )
 400      continue
C 
          incr = array(arr_p)*quotient    
          print *,  (array(entry + incr + j),
     +                    j = 0, remainder -1 )
      endif
      return
      end
C
C
C ******************************************************  
C  Reads values into local_entries member on each process.
C      If values should go into entries member, it
C      is necessary to add a loop to copy the values.
 
      subroutine Read_entries(prompt, array)
      character *21  prompt
      integer        array(0: 8 + 1024 + 1024-1)
C
      INCLUDE  'ciof.h'
      INCLUDE  'cyclic_iof.h'
      INCLUDE  'mpif.h'
      integer root
      integer i
      integer c
      integer recv_size
      integer retval
      integer ierr
C
      retval = Get_io_rank(array(arr_comm),  root)
C
      if (array(arr_myrank) .EQ. root)  then
          print *,'  ', prompt, '(one line per entry):'
          do 100 i = 0 ,array(global_order) - 1 
              read *,  array(entry + i)
 100      continue
          do 200 i = array(global_order), array(padded_size)-1,1
                  array(entry + i) = 0.0
                  print *,'pad updated!!'
 200      continue
C
      endif
C
      recv_size =  array(padded_size)/array(arr_p)
      call MPI_SCATTER(array(entry), 1, array(type),
     +     array(local_entry), recv_size, MPI_INTEGER ,
     +     root, array(arr_comm), ierr)
C
      return
      end  
