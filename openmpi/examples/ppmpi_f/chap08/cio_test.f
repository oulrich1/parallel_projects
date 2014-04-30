C  cio_test.f -- program for testing the functions in cio.c
C 
C  Input: An int, a float, and a string, a couple of times
C 
C  Output: Information on the communicators and input.
C
C  Note: Link with cio.o
C 
C  See Chap 8, pp. 142 & ff in PPMPI
C
      PROGRAM TstCIO
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
      integer io_comm_world
      integer even_comm
      integer odd_comm
      integer duped_comm
      integer split_key
      integer p
      integer my_rank
      integer ival
      real    fval
      integer sval
      integer ret_val
      integer ierr
      integer x
      character*15 source_comm_name
      character*15 new_comm_name
      character *40 prompt, title
C
        IO_KEY = MPI_KEYVAL_INVALID
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD, 
     +                           my_rank, ierr )
C
      print *,'Process ',my_rank,' > IO_KEY = ',
     +                IO_KEY
      call MPI_COMM_DUP(MPI_COMM_WORLD, 
     +                io_comm_world, ierr )
      x = Cache_io_rank(MPI_COMM_WORLD, 
     +                io_comm_world )
C
      print *,'Process ',my_rank,' > IO_KEY = ',  
     +                IO_KEY 
      source_comm_name = 'MPI_COMM_WORLD '
      new_comm_name = 'io_comm_world  '
      call Print_attr(source_comm_name, new_comm_name,
     +       io_comm_world)
C
      split_key = MOD(my_rank ,2 )
      if (split_key .EQ. 0)  then
          call MPI_COMM_SPLIT(io_comm_world, 
     +        split_key, my_rank,  even_comm, ierr)
          x = Cache_io_rank(io_comm_world, even_comm)
          source_comm_name = 'io_comm_world  '
          new_comm_name = 'even_comm      '
          call Print_attr(source_comm_name, new_comm_name,
     +           io_comm_world)
          x = Cache_io_rank(MPI_COMM_WORLD, even_comm )
          source_comm_name = 'MPI_COMM_WORLD '
          new_comm_name =    'even_comm      '
          call Print_attr(source_comm_name, new_comm_name,
     +           io_comm_world)
          call MPI_COMM_DUP(even_comm,  
     +               duped_comm, ierr)
          source_comm_name = 'even_comm      '
          new_comm_name =    'duped_comm     '
          call Print_attr(source_comm_name, new_comm_name,
     +           io_comm_world)
      else  
          call MPI_COMM_SPLIT(io_comm_world, 
     +         split_key, my_rank,  odd_comm, ierr)
          x = Cache_io_rank(io_comm_world, odd_comm)
          source_comm_name = 'io_comm_world  '
          new_comm_name =    'odd_comm       '
          call Print_attr(source_comm_name, new_comm_name,
     +           io_comm_world)
          x = Cache_io_rank(MPI_COMM_WORLD, odd_comm )
          source_comm_name = 'MPI_COMM_WORLD '
          new_comm_name =    'odd_comm       '
          call Print_attr(source_comm_name, new_comm_name,
     +           io_comm_world)
          call MPI_COMM_DUP(odd_comm,  
     +                 duped_comm, ierr)
          source_comm_name = 'odd_comm       '
          new_comm_name =    'duped_comm     '
          call Print_attr(source_comm_name, new_comm_name,
     +           io_comm_world)
      endif
C
      write(6,100) my_rank
 100  format(' ','Process ',i2,' > Before first call to Cread')
      prompt = 'Enter int,float,int: '     
      ret_val = Cread(MPI_COMM_WORLD, 3, prompt, ival, fval,sval)
      write(6,101) my_rank
 101  format(' ','Process ',i2,' > Call Cprint with MPI_COMM_WORLD')
      title = 'My values for int, float int: '
      ret_val = Cprint(MPI_COMM_WORLD, 3, title, ival,  fval, sval)
C
      write(6,102) my_rank
 102  format(' ','Process ',i2,' > Before second call to Cread')
      ret_val = Cread(io_comm_world, 3, prompt, ival,fval, sval)
      write(6,103) my_rank
 103  format(' ','Process ',i2,' > Call Cprint with io_comm_world')
      ret_val = Cprint(io_comm_world, 3, title, ival, fval, sval)
C
C      if (split_key .NE. 0) then
C          ival = 2*ival  
C          fval = 2*fval
C          write(6,104) my_rank
C 104      format(' ','Process ',i2,' > Call Cprint with odd_comm')
C          ret_val = Cprint(odd_comm, 3, title, ival,  fval,  sval)
C          ival = 2*ival
C          fval = 2*fval
C          write(6,105) my_rank
C 105      format(' ','Process ',i2,' > Call Cprint with duped_comm')
C          ret_val = Cprint(duped_comm, 3, title, ival,  fval, sval)
C      endif
C
      call MPI_FINALIZE(ierr)
      end
C
      subroutine Print_attr(source_comm_name, new_comm_name, 
     +                                    new_comm)
      character*(*)     source_comm_name
      character*(*)     new_comm_name
      integer            new_comm
C
      INCLUDE  'mpif.h'
      INCLUDE  'ciof.h'
      integer   io_rank
      integer   flag
      integer   my_rank
      integer   ierr
C
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      call MPI_ATTR_GET(new_comm, IO_KEY, io_rank, flag,ierr)
      if (flag .EQ. 0)  then
          write(6,100) my_rank, new_comm_name, source_comm_name
 100      format(' ','Process ',i2,
     +        ' > No attribute associated to IO_KEY in ',a15,
     +        ' from ',a15)
      else  
          write(6,200) my_rank, io_rank, new_comm_name, 
     +                 source_comm_name
 200      format(' ','Process ',i2,' > io_rank = ',i2,
     +        ' in ',a15,' (received from ',a15,')')
      endif
      return
      end

