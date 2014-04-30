C  cfopen.f -- open a file and write ints received from each process to it
C 
C  Input: none
C  Output: process ranks in file "testfile"
C 
C  Note: Link with cio.o
C 
C  See Chap 8, p. 157 in PPMPI
C 
      PROGRAM OpnFil
      INCLUDE 'mpif.h'
      INCLUDE 'ciof.h'
      integer       p
      integer       my_rank
      integer       io_rank
      integer       list(32)
      integer       io_comm
      integer       i
      integer       fp
      integer       MAX
      common        MAX
      integer       ierr
C
      integer Cfopen
C
      call MPI_INIT( ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD,  my_rank, ierr )
C
      MAX = 32
      call MPI_COMM_DUP(MPI_COMM_WORLD,  io_comm, ierr )
      call Cache_io_rank(MPI_COMM_WORLD, io_comm )
C
      fp = Cfopen('testfile',  io_comm)
C
      call Get_io_rank(io_comm,  io_rank)
      call MPI_GATHER( my_rank, 1, MPI_INTEGER, 
     +     list, 1, MPI_INTEGER,
     +     io_rank, io_comm, ierr)
      if (my_rank .EQ. io_rank)  then
          do 100 i = 1 , p
               write (unit = 3, FMT = 200) list(i)
 200           format (I7)
 100      continue
          end file (unit = 3)
          close (unit = 3)
      endif
C
      call MPI_FINALIZE(ierr)
      end
C
C
C **********************************************************  
      integer function Cfopen(filename, io_comm)
      character *25     filename     
      integer           io_comm    
C
      INCLUDE    'ciof.h'
      INCLUDE    'mpif.h'
      integer    root
      integer    my_io_rank
      integer    ierr
C
      call Get_io_rank(io_comm,  root)
C
      call MPI_COMM_RANK(io_comm,  my_io_rank, ierr)
C
      if (my_io_rank .EQ. root)  then
          OPEN (unit = 3, file = filename, form = 'FORMATTED' ,
     +                access = 'SEQUENTIAL')
          Cfopen = 1
      else  
          Cfopen = 0  
      endif
      return
      end   
