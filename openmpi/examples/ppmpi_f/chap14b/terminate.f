C  terminate.f -- functions for determining whether every process has
C      exhausted its stack.
C 
C  Here's the basic idea.  At the start of the program, process 0 has
C  one unit of some indestructible quantity -- I called it energy.
C  When the data is distributed, the energy is divided into p
C  equal parts and also distributed among the processes.  So each
C  process has 1/p of the energy.  (I don't actually send the energy
C  out initially, since every process knows it will get 1/p units of
C  energy, each can just set its energy to 1/p).  Whenever
C  a process exhausts its local stack, it sends whatever energy
C  it has to process 0.  Whenever a process satisfies a request
C  for work, it splits its energy in two equal pieces, keeping half
C  for itself and sending half to the process receiving the work.
C  Each process keeps track of its available energy in "my_energy".
C  Process 0, also keeps track of returned energy in "returned_energy."
C  Since energy is never destroyed, the sum of all energy on all the
C  processes will always be 1.  When process 0 finds that
C  "returned_energy" is 1, no process (including itself) will have
C  any work left, and it can broadcast a termination message to all
C  processes.  The only catch here is that floating point
C  arithmetic isn't exact.  So there are functions for
C  rational arithmetic:  a fraction a/b is represented as a pair of
C  longs (a,b).  Dividing by 2 changes (a,b) to (a,2b).  Adding
C  (a,b) + (c,d) = (ad + bc,bd).  This can still cause problems,
C  but there are no security checks.  Future work, no doubt.
C 
C  See Chap 14, pp. 334 & ff, in PPMPI.
C
C *******************************************************************  
      subroutine Print_divisors(p, my_rank)
      integer p
      integer my_rank
      integer i
      include 'terminatef.h'
C
      print 10, 'Process ',my_rank, ' > Divisors of ',2*p, 
     +          ' are : '
 10   format(A, I3, A, I3)
      print 20, ( divs(divisor_list +i), i = 0, divs(num_divisors)-1)
 20   format(I4)
      return
      end
C
C
C *******************************************************************  
      subroutine Setup_term_detect(p, my_rank)
      integer p
      integer my_rank
C
      include 'terminatef.h'
C
      call Find_divisors(2*p, p, my_rank)
C   Use 2*p to force inclusion of 2 in list   
C   of divisors                           
      my_energy(numerator) = 1
      my_energy(denominator) = p
C
      call Build_rational_mpi_t()
      return
      end   
C
C
C *******************************************************************  
      subroutine Build_rational_mpi_t()
      include 'terminatef.h'
      include 'mpif.h'
C
      integer       count
      integer       block_lengths(0:1)
      integer       types(0:1)
      integer       displacements(0:1)
      integer       start
      integer       address
      integer       ierr
C
      count = 2
      block_lengths(0) = 1
      block_lengths(1) = 1
      types(0) = MPI_INTEGER
      types(1) = MPI_INTEGER
C
      call MPI_ADDRESS( my_energy,  start, ierr)
      call MPI_ADDRESS( my_energy(numerator),  address, ierr)
      displacements(0) = address - start
      call MPI_ADDRESS( my_energy(denominator),  address, ierr)
      displacements(1) = address - start
C
      call MPI_TYPE_STRUCT(count, block_lengths, displacements, 
     +      types, rational_mpi_t, ierr)
      call MPI_TYPE_COMMIT( rational_mpi_t, ierr)
      return
      end    
C
C
C *******************************************************************  
      subroutine Find_divisors(  x, p, my_rank)
      integer x
      integer p
      integer my_rank
      include 'terminatef.h'  
      integer quotient
      integer remainder
      integer prime_indx
      integer divisor_indx
      integer   prime_list(0:MAX_PRIMES-1)
      data      prime_list
     +           /2,3,5,7,11,13,17,19,23,29,31,37,41,
     +            43,47,53,59,61,67,71,73,79,83,89,97/
      integer prime
C
      prime_indx= 0
      divisor_indx = 0
      quotient = x
      do while(quotient .NE. 1) 
          prime = prime_list(prime_indx)
          remainder = MOD(quotient , prime)
          if (remainder .EQ. 0)  then
              divs(divisor_list + divisor_indx) = prime
              divisor_indx = divisor_indx + 1
              quotient = quotient/prime
              remainder = MOD( quotient, prime)
              do while (remainder .EQ. 0)
                  quotient = quotient/prime
                  remainder = MOD( quotient, prime)
              end do
          end if
          prime_indx = prime_indx + 1
          if ((quotient .NE. 1)  .AND.  
     +        (prime_indx .GE. MAX_PRIMES))  then
              call Print_divisors(p, my_rank)
              print *,  'x = ',x, ' has too many distinct divisors '
              print *, 'Increase the size of prime_list' 
          endif
      end do
C
      divs(num_divisors) = divisor_indx
      return
      end
C
C
C *******************************************************************  
      subroutine Reduce( y) 
      integer y(0:1)
      include 'terminatef.h'
C
      integer   i
      integer   divisor
C
      do 100 i = 0 , divs(Num_divisors) -1   
          divisor =  divs(divisor_list +i)
          do while (MOD(y(numerator) , divisor) .EQ. 0
     +             .AND.
     +              MOD(y(denominator) , divisor) .EQ. 0) 
              y(Numerator) = y(Numerator)/divisor
              y(Denominator) = y(Denominator)/divisor
          end do
 100  continue
      return 
      end
C
C
C *******************************************************************  
C  Add x = x + y   
      subroutine Add( x,   y)  
      integer x(0:1)
      integer y(0:1)
      include 'terminatef.h'
C
      x(Numerator) = y(Denominator)*x(Numerator) + 
     +                x(Denominator)*y(Numerator)
      x(Denominator) = x(Denominator)*y(Denominator)
      call Reduce(x)
      return
      end
C
C
C *******************************************************************  
      subroutine Divide_by_2(x)
      integer x(0:1)
      include 'terminatef.h'   
C
      if (MOD(x(Numerator) , 2) .EQ. 0) then
          x(Numerator) = x(Numerator)/2
      else
          x(Denominator) = 2*x(Denominator)
      endif
      return
      end   
C
C
C *******************************************************************  
      integer function Equal( x, y) 
      integer x(0:1)
      integer y(0:1)
      include 'terminatef.h'
      include 'mainf.h'
C
      if ( x(Numerator)*y(Denominator) .EQ.
     +          y(Numerator)*x(Denominator) ) then
          Equal = TRUE
      else
          Equal = FALSE
      endif
      return
      end
C
C
C *******************************************************************  
      subroutine Send_half_energy(destination, comm)
      integer       destination   
      integer       comm , ierr         
      include 'terminatef.h'
      include 'mpif.h'
      include 'mainf.h'
C
      call Divide_by_2( my_energy)
C
      if (TREE_DEBUG) then
         print *,'Send Half Energy'
      endif
      call MPI_SEND( my_energy, 1, rational_mpi_t, destination,
     +     HALF_ENERGY_TAG, comm, ierr)
      return
      end   
C
C
C *******************************************************************  
      subroutine Recv_half_energy(source,comm)
      integer       source   
      integer       comm
      include       'terminatef.h'
      include       'mpif.h'     
      integer       status(MPI_STATUS_SIZE)
      integer       ierr
C
      call MPI_RECV( my_energy, 1, rational_mpi_t, source,
     +      HALF_ENERGY_TAG, comm,  status, ierr)
      return
      end
C
C
C *******************************************************************  
      subroutine Return_energy(comm) 
      integer comm
      integer ierr
      include 'terminatef.h'
      include 'mpif.h'
C
      call MPI_SEND( my_energy, 1, rational_mpi_t, 0, 
     +          RETURN_ENERGY_TAG, comm, ierr)
C
      return
      end   
C *******************************************************************  
      subroutine Receive_returned_energy(comm)
      integer         comm
      include 'terminatef.h'
      include 'mainf.h'   
      include 'mpif.h'
      integer         done
      integer         source
      logical         energy_returned
      integer         energy(0:1)
      integer         status(MPI_STATUS_SIZE)
C
      done = FALSE
      do while (done . EQ. FALSE) 
          call MPI_IPROBE(MPI_ANY_SOURCE, RETURN_ENERGY_TAG, 
     +        comm, energy_returned,  status, ierr)
          if (energy_returned )  then
              source = status(MPI_SOURCE)
              call MPI_RECV( energy, 1, rational_mpi_t, source,
     +             RETURN_ENERGY_TAG, comm,  status, ierr)
      if (TREE_DEBUG) then
          print 10,'Energy received = ',energy(0), energy(1)
          print 10, 'Add energy to = ', returned_energy(0),
     +                returned_energy(1)
 10       format(A,I3,I3)
      endif
              call Add( returned_energy, energy)
          else  
              done = TRUE
          endif
      end do
      return
      end  
C
C
C *******************************************************************  
      integer function Search_complete(comm, p, my_rank)
      integer     comm   
      integer     dest
      integer     x
      logical     completed
      include  'terminatef.h'
      include  'mainf.h'
      include  'mpif.h'
      integer     status(MPI_STATUS_SIZE)
C
      x = 0
      if (my_rank .EQ. 0)  then
          call Receive_returned_energy(comm)
          if (Equal( returned_energy,  ONE) .EQ. TRUE)  then
              do 100 dest = 1, p-1
                  call MPI_SEND( x, 1, MPI_INTEGER, dest, 
     +                 COMPLETE_TAG, comm, ierr)
 100          continue
              Search_complete =  TRUE
          else  
              Search_complete =  FALSE
          endif
      else
          call MPI_IPROBE(0, COMPLETE_TAG, comm,
     +                   completed,  status, ierr)
          if (completed  )  then
              call MPI_RECV(  x, 1, MPI_INTEGER, 0, COMPLETE_TAG, 
     +               comm,  status, ierr)
              Search_complete = TRUE
          else  
              Search_complete = FALSE
          endif
      endif
      if (TREE_DEBUG .AND. Search_complete .EQ. TRUE) then
          print *,my_rank, ' Search Complete!'
      endif
C
      return
      end   
