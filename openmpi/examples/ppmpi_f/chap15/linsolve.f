C   linsolve.f
C   Use Scalapack and MPI to solve a system of linear equations
C   on a virtual rectangular grid of processes.
C
C   Input:
C       N: order of linear system
C       NPROC_ROWS: number of rows in process grid
C       NPROC_COLS: number of columns in process grid
C       ROW_BLOCK_SIZE:  blocking size for matrix rows
C       COL_BLOCK_SIZE:  blocking size for matrix columns
C
C   Output:
C       Input data, error in solution, and time to solve system.
C
C   Algorithm:
C	1.  Initialize MPI and BLACS.
C       2.  Get process rank (MY_RANK) and total number of 
C           processes (NP).   Use both BLACS and MPI.
C       3a. Process 0 read and broadcast matrix order (N),
C           number of process rows (NPROC_ROWS), number 
C           of process columns (NPROC_COLS), ROW_BLOCK_SIZE,
C	    and COL_BLOCK_SIZE.
C       3b. Process != 0 receive same.
C       4.  Use BLACS_GRIDINIT to set up process grid.
C       5.  Compute amount of storage needed for local arrays.
C       6.  Use ScaLAPACK routine DESCINIT to initialize
C           descriptors for A, EXACT (= exact solution), and B.
C       7.  Use random number generator to generate contents 
C           of local block of matrix (A_LOCAL).
C       8.  Set entries of EXACT to 1.0.
C	9.  Generate B by computing B = A*EXACT.  Use
C           PBLAS routine PSGEMV.
C      10.  Solve linear system by call to ScaLAPACK routine
C           PSGESV (solution returned in B). 
C      11.  Use PBLAS routines PSAXPY and PSNRM2 to compute
C           the norm of the error ||B - EXACT||_2.
C      12.  Process 0 print results.
C      13.  Free up storage, shutdown BLACS and MPI.
C
C   Notes:
C       1.  The vectors EXACT, and B are significant only
C           in the first process column.
C       2.  A_LOCAL is allocated as a linear array.
C       3.  The solver only allows square blocks.  So we
C           read in a single value, BLOCK_SIZE, and assign
C           it to ROW_BLOCK_SIZE and COL_BLOCK_SIZE.  Thus,
C           since the matrix is square, NPROC_ROWS = NPROC_COLS.
C
	PROGRAM LINSOLVE
	INCLUDE 'mpif.h'
C
C   Constants
	INTEGER		MAX_VECTOR_SIZE
	INTEGER		MAX_MATRIX_SIZE
	INTEGER		DESCRIPTOR_SIZE
        PARAMETER	(MAX_VECTOR_SIZE = 1000)
	PARAMETER	(MAX_MATRIX_SIZE = 250000)
	PARAMETER	(DESCRIPTOR_SIZE = 8)
C
C   Array Variables
	REAL		B_LOCAL(MAX_VECTOR_SIZE)
	INTEGER		B_DESCRIP(DESCRIPTOR_SIZE)
C
	REAL		EXACT_LOCAL(MAX_VECTOR_SIZE)
	INTEGER		EXACT_DESCRIP(DESCRIPTOR_SIZE)
C
	REAL 		A_LOCAL(MAX_MATRIX_SIZE)
	INTEGER		A_DESCRIP(DESCRIPTOR_SIZE)
	INTEGER		PIVOT_LIST(MAX_VECTOR_SIZE)
C
C   Scalar Variables
	INTEGER 	NP
	INTEGER		MY_RANK
	INTEGER		NPROC_ROWS
	INTEGER		NPROC_COLS
	INTEGER		IERROR
	INTEGER		M, N
	INTEGER		ROW_BLOCK_SIZE
	INTEGER		COL_BLOCK_SIZE
        INTEGER		INPUT_DATA_TYPE
        INTEGER		BLACS_CONTEXT
        INTEGER		TEMP_CONTEXT(1)
        INTEGER		LOCAL_MAT_ROWS
        INTEGER		LOCAL_MAT_COLS
	INTEGER		EXACT_LOCAL_SIZE
	INTEGER		B_LOCAL_SIZE
	INTEGER		I, J
	INTEGER		MY_PROCESS_ROW
	INTEGER		MY_PROCESS_COL
        REAL		ERROR_2
	DOUBLE PRECISION START_TIME
	DOUBLE PRECISION ELAPSED_TIME
C
C   Local Functions
C       REAL		RAND_VAL
C
C   External subroutines
C     MPI:
	EXTERNAL	MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK
	EXTERNAL	MPI_BCAST, MPI_ABORT
	EXTERNAL	MPI_BARRIER, MPI_FINALIZE
C     BLACS:
	EXTERNAL	BLACS_GET, BLACS_GRIDINIT, BLACS_PCOORD
	EXTERNAL	BLACS_EXIT
C     PBLAS:
	EXTERNAL	PSGEMV, PSAXPY, PSNRM2
C     ScaLAPACK
	EXTERNAL	DESCINIT, PSGESV
C     System
	EXTERNAL	SLARNV
C
C   External functions
C     ScaLAPACK:
	INTEGER		NUMROC
	EXTERNAL	NUMROC
C
C   Junk Variables
	INTEGER		ISEED(4)
C
C   Begin Executable Statements
C
C   Initialize MPI and BLACS
	CALL MPI_INIT(IERROR)
C
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, IERROR)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, IERROR)
C	
C
C   Get Input Data.  First Build Derived Type.
	CALL BUILD_INPUT_DATA_TYPE(INPUT_DATA_TYPE, N, NPROC_ROWS,
     +                             NPROC_COLS, ROW_BLOCK_SIZE,
     +                             COL_BLOCK_SIZE)
	IF (MY_RANK.EQ.0) THEN
	    READ(5,*) N, NPROC_ROWS, NPROC_COLS, ROW_BLOCK_SIZE,
     +                COL_BLOCK_SIZE
	END IF
	CALL MPI_BCAST(N, 1, INPUT_DATA_TYPE, 0, MPI_COMM_WORLD, IERROR)
	IF (NP.LT.(NPROC_ROWS*NPROC_COLS)) THEN
	    WRITE(6,250) MY_RANK, NP, NPROC_ROWS, NPROC_COLS
  250	    FORMAT(' ','Proc ',I2,' > NP = ',I2,', NPROC_ROWS = ',I2,
     +             ', NPROC_COLS = ',I2)
	    WRITE(6,260)
  260       FORMAT(' ','Need more processes!  Quitting.')
	    CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
        END IF
C
C   The matrix is square
        M = N
C
C
C   Build BLACS grid.
C     First get BLACS System Context (in TEMP_CONTEXT(1))
	CALL BLACS_GET(0, 0, TEMP_CONTEXT)
	BLACS_CONTEXT = TEMP_CONTEXT(1)
C
C     BLACS_CONTEXT is in/out.
C     'R': process grid will use row major ordering.
	CALL BLACS_GRIDINIT(BLACS_CONTEXT,'R',NPROC_ROWS,
     +                        NPROC_COLS)
C
C
C   Figure out how many rows and cols we'll need in the local
C     matrix.
	CALL BLACS_PCOORD(BLACS_CONTEXT, MY_RANK, MY_PROCESS_ROW,
     +                    MY_PROCESS_COL)
	LOCAL_MAT_ROWS = NUMROC(M, ROW_BLOCK_SIZE, MY_PROCESS_ROW,
     +				0, NPROC_ROWS)
	LOCAL_MAT_COLS = NUMROC(N, COL_BLOCK_SIZE, MY_PROCESS_COL,
     +			        0, NPROC_COLS)
	IF (LOCAL_MAT_ROWS*LOCAL_MAT_COLS.GT.MAX_MATRIX_SIZE) THEN
	    WRITE(6,290) MY_RANK, LOCAL_MAT_ROWS, LOCAL_MAT_COLS,
     +                   MAX_MATRIX_SIZE
  290       FORMAT(' ','Proc ',I2,' > LOCAL_MAT_ROWS = ',I5,
     +             ', LOCAL_MAT_COLS = ',I5,
     +             ', MAX_MATRIX_SIZE = ',I6)
	    WRITE(6,292)
  292       FORMAT(' ','Insufficient storage!  Quitting.')
	    CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
	END IF
C
C   Now figure out storage for B_LOCAL and EXACT_LOCAL
	B_LOCAL_SIZE = NUMROC(M, ROW_BLOCK_SIZE, MY_PROCESS_ROW,
     +			      0, NPROC_ROWS)
	IF (B_LOCAL_SIZE.GT.MAX_VECTOR_SIZE) THEN
	    WRITE(6,294) MY_RANK, B_LOCAL_SIZE, MAX_VECTOR_SIZE
  294       FORMAT(' ','Proc ',I2,' > B_LOCAL_SIZE = ',I5,
     +             ', MAX_VECTOR_SIZE = ',I5)
	    WRITE(6,296)
  296       FORMAT(' ','Insufficient storage!  Quitting.')
	    CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
	END IF
C
	EXACT_LOCAL_SIZE = NUMROC(N, COL_BLOCK_SIZE, MY_PROCESS_ROW,
     +			      0, NPROC_ROWS)
	IF (EXACT_LOCAL_SIZE.GT.MAX_VECTOR_SIZE) THEN
	    WRITE(6,298) MY_RANK, EXACT_LOCAL_SIZE, MAX_VECTOR_SIZE
  298       FORMAT(' ','Proc ',I2,' > EXACT_LOCAL_SIZE = ',I5,
     +             ', MAX_VECTOR_SIZE = ',I5)
	    WRITE(6,299)
  299       FORMAT(' ','Insufficient storage!  Quitting.')
	    CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
	END IF
C
C
C   Now build the matrix descriptors using ScaLAPACK'S DESCINIT
	CALL DESCINIT(A_DESCRIP, M, N, ROW_BLOCK_SIZE, 
     +                COL_BLOCK_SIZE, 0, 0, BLACS_CONTEXT,
     +                LOCAL_MAT_ROWS, IERROR)
	IF (IERROR.NE.0) THEN
            WRITE(6,300) MY_RANK, IERROR 
  300       FORMAT(' ','Proc ',I2,
     +             ' > DESCINIT FOR A FAILED, IERROR = ',I3)
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
        END IF
C
        CALL DESCINIT(B_DESCRIP, M, 1, ROW_BLOCK_SIZE, 1,
     +                0, 0, BLACS_CONTEXT, B_LOCAL_SIZE,
     +                IERROR)
	IF (IERROR.NE.0) THEN
            WRITE(6,350) MY_RANK, IERROR 
  350       FORMAT(' ','Proc ',I2,
     +		   ' > DESCINIT FOR B FAILED, IERROR = ',I3)
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
        END IF
C       
        CALL DESCINIT(EXACT_DESCRIP, N, 1, COL_BLOCK_SIZE, 1,
     +                0, 0, BLACS_CONTEXT, EXACT_LOCAL_SIZE,
     +                IERROR)
	IF (IERROR.NE.0) THEN
            WRITE(6,400) MY_RANK, IERROR 
  400       FORMAT(' ','Proc > ',I2,
     +		   ' > DESCINIT FOR EXACT FAILED, IERROR = ',I3)
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
        END IF
C
C
C   Now initialize A_LOCAL and EXACT_LOCAL
C	CALL SEED(MY_RANK)
C	DO 600 J = 0, LOCAL_MAT_COLS - 1
C	    DO 500 I = 1, LOCAL_MAT_ROWS
C		A_LOCAL(LOCAL_MAT_ROWS*J + I) = RAND_VAL()
C  500	    CONTINUE
C  600  CONTINUE
C     SGI subroutine -- generates the full matrix
	ISEED(1) = MY_RANK
        ISEED(2) = MY_RANK*MY_RANK
        ISEED(3) = NP - MY_RANK
        ISEED(4) = 2*MY_RANK + 1
	DO 650 J = 0, LOCAL_MAT_COLS - 1
            CALL SLARNV(1, ISEED, LOCAL_MAT_ROWS, 
     +        		A_LOCAL(J*LOCAL_MAT_ROWS + 1))
  650	CONTINUE
C
	DO 700 I = 1, EXACT_LOCAL_SIZE
	    EXACT_LOCAL(I) = 1.0
  700	CONTINUE
C
C
C   Use PBLAS function PSGEMV to compute right-hand side B = A*EXACT
C     'N': Multiply by A -- not A^T or A^H
	CALL PSGEMV('N', M, N, 1.0, A_LOCAL, 1, 1, A_DESCRIP,
     +              EXACT_LOCAL, 1, 1, EXACT_DESCRIP, 1, 0.0,
     +              B_LOCAL, 1, 1, B_DESCRIP, 1)
C 
C
C   Done with setup!  Solve the system.
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
	START_TIME = MPI_WTIME()
	CALL PSGESV(N, 1, A_LOCAL, 1, 1, A_DESCRIP, PIVOT_LIST,
     +              B_LOCAL, 1, 1, B_DESCRIP, IERROR)
	ELAPSED_TIME = MPI_WTIME() - START_TIME
	IF (IERROR.NE.0) THEN
            WRITE(6,800) MY_RANK, IERROR 
  800       FORMAT(' ','Proc ',I2,' > PSGESV FAILED, IERROR = ',I3)
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, IERROR)
        END IF
C
C	WRITE(6,850) MY_RANK, (B_LOCAL(J), J = 1, B_LOCAL_SIZE)
C 850   FORMAT(' ','Proc ',I2,' > B = ',F6.3,' ',F6.3,' ',F6.3,' ',
C    +         F6.3,' ',F6.3,' ',F6.3,' ',F6.3,' ',F6.3,' ',F6.3,' ',
C    +         F6.3)
C	WRITE(6,860) MY_RANK, (EXACT_LOCAL(I), I = 1, EXACT_LOCAL_SIZE)
C 860   FORMAT(' ','Proc ',I2,' > EXACT = ',F6.3,' ',F6.3,' ',F6.3,' ',
C    +         F6.3,' ',F6.3,' ',F6.3,' ',F6.3,' ',F6.3,' ',F6.3,' ',
C    +         F6.3)
C
C   Now find the norm of the error.
C     First compute EXACT = -1*B + EXACT
        CALL PSAXPY(N, -1.0, B_LOCAL, 1, 1, B_DESCRIP, 1,
     +              EXACT_LOCAL, 1, 1, EXACT_DESCRIP, 1)
C     Now compute 2-norm of EXACT
	CALL PSNRM2(N, ERROR_2, EXACT_LOCAL, 1, 1, EXACT_DESCRIP, 1)
C
C
	IF (MY_RANK.EQ.0) THEN
	    WRITE(6,900) N, NP
  900       FORMAT(' ','N = ',I4,', Number of Processes = ',I2)
	    WRITE(6,950) NPROC_ROWS, NPROC_COLS
  950	    FORMAT(' ','Process rows = ',I2,', Process cols = ',I2)
	    WRITE(6,1000) ROW_BLOCK_SIZE, COL_BLOCK_SIZE
 1000	    FORMAT(' ','Row block size = ',I3,', Col block size = ',I3)
	    WRITE(6,1100) ERROR_2
 1100	    FORMAT(' ','2-Norm of error = ',E13.6)
	    WRITE(6,1200) 1000.0*ELAPSED_TIME
 1200	    FORMAT(' ','Elapsed time = ',D13.6,' milliseconds')
        END IF
C
C
C   Now free up allocated resources and shut down
C     Call BLACS_EXIT.  Argument != 0 says, "I'll shut down MPI."
	CALL BLACS_EXIT(1)
        CALL MPI_FINALIZE(IERROR)
C
	STOP
C
C 	End of Main Program LINSOLVE
	END
C
C
C**********************************************************************
 	SUBROUTINE BUILD_INPUT_DATA_TYPE(INPUT_DATA_TYPE, N, 
     +                                   NPROC_ROWS, NPROC_COLS, 
     +                                   ROW_BLOCK_SIZE, 
     +                                   COL_BLOCK_SIZE)
C 
C  Build a derived datatype for transmitting the five input 
C  values in a single broadcast 
C 
	INCLUDE 'mpif.h'
C
C  Local Constant
	INTEGER  	ELEMENTS 
	PARAMETER 	(ELEMENTS = 5) 
C 
C  Output Parameter 
	INTEGER 	INPUT_DATA_TYPE 
C 
C  Input Parameters 
	INTEGER 	N
	INTEGER 	NPROC_ROWS 
	INTEGER  	NPROC_COLS 
	INTEGER  	ROW_BLOCK_SIZE 
	INTEGER  	COL_BLOCK_SIZE 
C
C  Local Array Variables 
	INTEGER 	ARRAY_OF_BLOCK_LENGTHS(ELEMENTS)
	INTEGER 	ARRAY_OF_DISPLACEMENTS(ELEMENTS) 
	INTEGER		ARRAY_OF_TYPES(ELEMENTS) 
C
C  Local Scalar Variables 
	INTEGER 	BASE_ADDRESS 
	INTEGER  	TEMP_ADDRESS 
	INTEGER  	IERROR 
	INTEGER		I
C
C   External Subroutines
C     MPI:
	EXTERNAL	MPI_ADDRESS, MPI_TYPE_STRUCT, MPI_TYPE_COMMIT
C 
	DATA ARRAY_OF_BLOCK_LENGTHS / 1, 1, 1, 1, 1/ 
C 
	DO 100 I = 1, ELEMENTS 
     	    ARRAY_OF_TYPES(I) = MPI_INTEGER 
  100 	CONTINUE 
C 
C   Compute displacements from N 
	ARRAY_OF_DISPLACEMENTS(1) = 0 
	CALL MPI_ADDRESS(N, BASE_ADDRESS, IERROR) 
	CALL MPI_ADDRESS(NPROC_ROWS, TEMP_ADDRESS, IERROR) 
	ARRAY_OF_DISPLACEMENTS(2) = TEMP_ADDRESS - BASE_ADDRESS
	CALL MPI_ADDRESS(NPROC_COLS, TEMP_ADDRESS, IERROR)
	ARRAY_OF_DISPLACEMENTS(3) = TEMP_ADDRESS - BASE_ADDRESS
	CALL MPI_ADDRESS(ROW_BLOCK_SIZE, TEMP_ADDRESS, IERROR)
	ARRAY_OF_DISPLACEMENTS(4) = TEMP_ADDRESS - BASE_ADDRESS 
	CALL MPI_ADDRESS(COL_BLOCK_SIZE, TEMP_ADDRESS, IERROR)
	ARRAY_OF_DISPLACEMENTS(5) = TEMP_ADDRESS - BASE_ADDRESS 
C
	CALL MPI_TYPE_STRUCT(ELEMENTS, ARRAY_OF_BLOCK_LENGTHS, 
     +                  ARRAY_OF_DISPLACEMENTS, 
     +                  ARRAY_OF_TYPES, 
     +                  INPUT_DATA_TYPE, IERROR) 
	CALL MPI_TYPE_COMMIT(INPUT_DATA_TYPE, IERROR) 
C
	RETURN 
C 
C   End of Subroutine BUILD_INPUT_DATA_TYPE
	END 
C
C
C**********************************************************************
