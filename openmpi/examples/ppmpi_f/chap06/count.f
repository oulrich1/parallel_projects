C  count.f -- send a subvector from process 0 to process 1
C 
C  Input: none
C  Output: contents of vector received by process 1
C 
C  Note: Program should only be run with 2 processes.
C 
C  See Chap 6, pp. 89 & ff. in PPMPI
C  
C
        program count
        include 'mpif.h'
C
        real vector(100)
        integer status(MPI_STATUS_SIZE)
        integer p
        integer my_rank
        integer i
        integer ierr

        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr)
    	call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

C  Initialize vector and send */
    	if (my_rank .EQ. 0) then
        	do 100 i = 1, 50
            		vector(i) = 0.0
 100		continue
        	do 200 i = 51, 100  
          		vector(i) = 1.0
 200		continue
        	call MPI_SEND(vector(51), 50, MPI_REAL, 1, 0,
     +       		MPI_COMM_WORLD, ierr) 
    	else  
        	call MPI_RECV(vector(51), 50, MPI_REAL, 0, 0,
     +                          MPI_COMM_WORLD, status, ierr)
        	do 300 i = 51, 100
            		print *, vector(i), ' '
 300		continue
        	print *, ' '
    	endif

    	call MPI_FINALIZE(ierr)
        end
