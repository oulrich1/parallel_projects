C 
C  terminatef.h
C 
C  Definitions and declarations for distributed termination detection
C
        integer RETURN_ENERGY_TAG  
        integer HALF_ENERGY_TAG    
        integer COMPLETE_TAG 
        parameter (RETURN_ENERGY_TAG = 1000,
     +              HALF_ENERGY_TAG = 2000,
     +              COMPLETE_TAG = 3000)
C  Number of primes stored in prime_list   
C      Currently number of primes < 100    
        integer MAX_PRIMES 
        parameter (MAX_PRIMES = 25)
C
C  typedef struct  DIVISORS
      integer   num_divisors
      integer   divisor_list
      integer   sizeofDIVS
      parameter (num_divisors = 0, divisor_list = 1,
     +            sizeofDIVS = 26)
      integer   divs
      common /Term/   divs(0:sizeofDIVS-1)
C
C  typedef struct RATIONAL
      integer rational_mpi_t
      common /RatlMPI/ rational_mpi_t
      integer numerator
      integer denominator
      parameter (numerator = 0, denominator = 1)
      integer  ONE(0:1) 
      integer  my_energy(0:1)
      integer  returned_energy(0:1)
      common /Enrgy/ my_energy, returned_energy,
     +                ONE
C
C
C  functions
      integer   Equal 
      integer   Search_complete 
C
