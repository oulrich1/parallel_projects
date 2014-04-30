 /* 
    Matrix Mult in Parallel:

    By Oriah Ulrich

    Accepts as input, two matricies and multiplies them together
    the order of the matrix must be N x N square matrix, since
    the program has not yet been optimized to handle edge cases. 
    With that said, N must be a multiple of the number of processors
    used to run the program. With p processors N must divide evenly 
    into p processors..

    Using premake to build the Makefiles, but the program can also
    be compiled with the following command. The following command
    takes advantage of gcc compiler optimizations and processor
    optimizations. 

    To find the working machines and generate the hostfile, just run
    ./find_machines in the directory where you would like to put the 
    hostfiles.. (Since there are issues with hanging runs on jaguar when
        trying to run on machines that cannot be connected to)


    To compile use mpicc:

            mpicc -g -Wall -std=c99  \
                    -O3 -funroll-loops -march=native  \
                    -o matrixmult ./matrixmult.c  \
                    -Wno-unused-result  \
                    -Wno-unused-variable  -lm

    To run, use mpirun:

        mpirun -hostfile 4hosts_2procs ./matrixmult < tests/t4800.out

        mpirun ./matrixmult < tests/t3000.ijk.out

    Design description
        At first the idea was to perform a parallel dot product
        but after looking at the order of convergence it seemed 
        like it would actually be slower by a factor of log p

        Instead i chose to parallelize the data sent to the processor
        However, in order to get this to work with any sized matrix and
        handle all edge cases i need to consider:


            //does N divide into the number of processors 
            if (N < proc_count) {
                // assign a processor to each row
                // calculate remaining processors
                rows_pproc = 1;
                rem_rows  = 0;
                rem_procs = proc_count - N;
                printf("Still need to handle extra processors...
                        proccessors left to work: %d\n",  proc_count);
            } else if (N == proc_count) {
                rows_pproc = 1;
                rem_rows  = 0;
                rem_procs = 0;
                // this works out nicely with number of procs, 
                // but this is not likely to happen anyways
            } else if (N > proc_count) {
                most likely case when N > proc_count:
                    will assign blocks of rows to processors
                    each block will be the same size of elements
                    
                rows_pproc = N / proc_count;
                rem_rows  = N % proc_count;
                rem_procs = 0;

                THIS IS VERY IMPORTANT!
                    this represents the id of the row (0 indexed) 
                    where we need to append a new partial matrix 
                    mult solution to. 
                unassigned_row_id = (rows_pproc * proc_count); 
            }  
                
        However part has not yet been fulling implemented, more work is
        necessary to ensure the data is passed back properly and dynamically..

        Fow now, the program makes the assumption that N will be a multiple of
        the number of cores being utilized..

        As the goal of the project is to look into parallelizing the matrix mult 
        for ijk forms on a 4800 x 4800 matrix, i will implement the remaining 
        partition of the column ordered transposed matrix later...

*/

    /* MUST COMPILE WITH -lm */

#include <time.h>

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdarg.h>
// #include <unistd.h> /* machine info */

#include <string.h>
#include <math.h>


#ifndef NULL
#define NULL 0
#endif

typedef unsigned int uint;
typedef long double LD_T;

typedef double PRIM_T;

#include <mpi.h>

#define MatrixEntry(A,i,j,ne) (*( (A) + (ne)*(i) + (j) ))

#define STDIN_PROMPT_ARGS       false

// set this to true if you want
// to see the final result
#define STDIN_RESULT_OUTPUT     false   


#define OUTPUT_HOSTNAME         true


enum matrix_mult_state_t
{
    NULL_MODE = 0,
    IJK_MODE  = 1,
    IKJ_MODE  = 2,
    KIJ_MODE  = 4
};

enum matrix_genr_state_t
{
    INPUT_MODE  = 1,
    RAND_MODE  = 2
};

#define MASTER 0

void matrix_print(double* A, int m, int n) {
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            printf("  %.2lf ", A[ i * n + j ]);
        }
        printf("\n");
    }
    printf("\n");
}


int main(int argc, char const *argv[])
{       
    int   comm_sz;               /* Number of processes    */
    int   my_rank;               /* My process rank        */

    /* timing variables */
    double time_local;
    double elapsed;

    int   N = 100;
    int   err = 0; (void) err;
    char  state_str[4] = "ijk";
    char  gen_flag[2] = "I";
    uint  matrix_mult_state = IJK_MODE; // ijk, ikj, kij forms
    uint  matrix_genr_state = RAND_MODE;
    
    double* A;
    double* B;       // B is local to all processors
    double* C;
    double* A_REM = NULL;
    double* C_REM = NULL;

    /* Start up MPI */
    MPI_Init(NULL, NULL); 

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 

    /* Get my rank among all the processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

    /* Verify we are on the right computers.. */
    #if OUTPUT_HOSTNAME == true
    for (int i = 0; i < comm_sz; ++i) {
        if (my_rank == i) {
            char hostname[1024];
            gethostname(hostname, 1024);
            printf("My Rank: #%d. Running on machine: %s\n", my_rank, hostname);
            MPI_Barrier(MPI_COMM_WORLD);  
        } else {
            MPI_Barrier(MPI_COMM_WORLD); 
        }
    }
    #endif

    if (MASTER == my_rank)
    {
        printf("------------------------\n");
        printf("running on %d processors\n", comm_sz);
    }
    
    /* Collect STDIN data:
        ijk form
        input flag
        Matrix Size N */
    if (MASTER == my_rank)
    {
        #if STDIN_PROMPT_ARGS == true
        printf("input multiply form: ");
        #endif
        scanf("%s", state_str);

        if (!strcmp(state_str, "ijk")) { 
            matrix_mult_state = IJK_MODE;
        } else if (!strcmp(state_str, "ikj")) { 
            matrix_mult_state = IKJ_MODE;
        } else if (!strcmp(state_str, "kij")) { 
            matrix_mult_state = KIJ_MODE;
        } else {
            printf("unknown matrix mode.. defaulting to ijk.\n");
            matrix_mult_state = IJK_MODE;
        }

        #if STDIN_PROMPT_ARGS == true
        printf("input flag: [ R - Random Matrix | I - Stdin ]:");
        #endif
        scanf("%s", gen_flag);
        if (!strcmp(state_str, "I") || !strcmp(state_str, "i")) { 
            matrix_genr_state = INPUT_MODE;
        } else if (!strcmp(state_str, "R") || !strcmp(state_str, "r")) { 
            matrix_genr_state = RAND_MODE;
        }

        #if STDIN_PROMPT_ARGS == true
        printf("input square matrix size: ");
        #endif
        scanf("%d", &N);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    time_local = MPI_Wtime();

    err = MPI_Bcast(&matrix_mult_state, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
    err = MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    time_local = MPI_Wtime() - time_local;

    int proc_count  = comm_sz;
    int rem_procs   = 0; 
        (void) rem_procs;
    int rows_pproc  = 0;
    int rem_rows    = 0;
    int unassigned_row_id = -1;
        (void) unassigned_row_id;

    /* does N divide into the number of processors */
    if (N < proc_count) {
        // assign a processor to each row
        // calculate remaining processors
        rows_pproc = 1;
        rem_rows  = 0;
        rem_procs = proc_count - N;
        printf("Still need to handle extra processors... proccessors left to work: %d\n",  proc_count);
    } else if (N == proc_count) {
        rows_pproc = 1;
        rem_rows  = 0;
        rem_procs = 0;
        // this works out nicely with number of procs, 
        // but this is not likely to happen anyways
    } else if (N > proc_count) {
        /* most likely case when N > proc_count:
            will assign blocks of rows to processors
            each block will be the same size of elements
             */
        rows_pproc = N / proc_count;
        rem_rows  = N % proc_count;
        rem_procs = 0;

        /* THIS COULD BE VERY IMPORTANT! (NOT YET FULLY IMPLEMENTED YET)
            this represents the id of the row (0 indexed) 
            where we need to append a new partial matrix 
            mult solution to. */
        unassigned_row_id = (rows_pproc * proc_count); 
    }

    /* parallelized and gathered matricies' matrix size variables */

    int A_row_count = rows_pproc; // each processor gets this many rows..
    int A_local_row_count = rows_pproc; // each processor gets this many rows..
        (void) A_local_row_count;

    int B_col_count = N;

    int A_col_count = N;
    int B_row_count = N;

    int vector_size = N;

    /* allocate the matricies */

    if (my_rank == MASTER) {
        A = (double*) malloc(sizeof(double) * N * N);
        C = (double*) malloc(sizeof(double) * N * N);
    } else {
        /* each proc gets the same number of elements..
            if the proc number  */

        int count_elements_in_A = N * rows_pproc;
        A = (double*) malloc(sizeof(double) * count_elements_in_A);
        C = (double*) malloc(sizeof(double) * count_elements_in_A); 
        
        /* these processors will also carry on the
            burden of calculating the residue matrix 
            multiplication */
        if (my_rank <= rem_rows-1) {
            A_REM = (double*) malloc(sizeof(double) * N);
            C_REM = (double*) malloc(sizeof(double) * N);
        } else {
            A_REM = NULL;
            C_REM = NULL;
        }
    }

    B = (double*) malloc(sizeof(double) * N * N);


    /* Populate the matricies with values from sources:
        STDIN or Randomly Generated */

    /* If GENR MODE IS INPUT MODE then :
            Collect STDIN data:
                ijk form
                input flag
                Matrix Size N */
    if (my_rank == MASTER)
    {
        double element_value = 0;
        /*const*/ int MAX_RAND_VALUE = 100;
        /*const*/ int POSITIVE_FACTOR = 2;

        /*      INITIALIZE THE MATRICIES WITH THE STD INPUT     */
        /* First if else block reads in Matrix A */
        /* Second if else block reads in Matrix B */

        /* this is the first if else block: Read in Matrix A */
        if (matrix_mult_state & KIJ_MODE)
        // if(false)
        {
             // set to this if the kij form is crashing.. 
             // and change the iteration of A in kij
            
        /* if KIJ then store the Matrix B as a transpose of itself.. 
            (for faster access and less cache miss) */
            #if STDIN_PROMPT_ARGS == true
            printf("\ninput %d values for matrix: #%d (row major form) \n", N * N, 1);
            #endif
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (matrix_genr_state & RAND_MODE) {
                        double random  = sqrt((rand() % MAX_RAND_VALUE)) 
                                                * (rand() % MAX_RAND_VALUE);
                        int neg_or_pos = ((2 * (rand() % POSITIVE_FACTOR)) - 1);
                        element_value = random * neg_or_pos;
                    } else if (matrix_genr_state & INPUT_MODE) {
                        scanf("%lf", &element_value);
                    }
                    A[ j * N + i ] = element_value;
                }
            }
        } else {
            #if STDIN_PROMPT_ARGS == true
            printf("\ninput %d values for matrix: #%d (row major form) \n", N * N, 1);
            #endif
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (matrix_genr_state & RAND_MODE) {
                        double random  = sqrt((rand() % MAX_RAND_VALUE)) 
                                                * (rand() % MAX_RAND_VALUE);
                        int neg_or_pos = ((2 * (rand() % POSITIVE_FACTOR)) - 1);
                        element_value = random * neg_or_pos;
                    } else if (matrix_genr_state & INPUT_MODE) {
                        scanf("%lf", &element_value);
                    }
                    A[ i * N + j ] = element_value;
                }
            }
        }

        /* this is the second if else block: Read in Matrix B */
        /* if IJK then store the Matrix B as a transpose of itself.. 
            (for faster access and less cache miss) */
        if (matrix_mult_state & IJK_MODE) {
            #if STDIN_PROMPT_ARGS == true
            printf("\ninput %d values for matrix: #%d (column major form) \n\n", N * N, 2);
            #endif
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (matrix_genr_state & RAND_MODE) {
                        double random  = sqrt((rand() % MAX_RAND_VALUE)) 
                                                * (rand() % MAX_RAND_VALUE);
                        int neg_or_pos = ((2 * (rand() % POSITIVE_FACTOR)) - 1);
                        element_value = random * neg_or_pos;
                    } else if (matrix_genr_state & INPUT_MODE) {
                        scanf("%lf", &element_value);
                    }
                    // B[ i * N + j ] = element_value;
                    B[ j * N + i ] = element_value; // make each row element a column element (transpose)
                }
            }
        } else {
            /* else it is ikj or kij, 
                can do row major storage of B */
            #if STDIN_PROMPT_ARGS == true
            printf("\ninput %d values for matrix: #%d (row major form) \n\n", N * N, 2);
            #endif
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (matrix_genr_state & RAND_MODE) {
                        double random  = sqrt((rand() % MAX_RAND_VALUE)) 
                                                * (rand() % MAX_RAND_VALUE);
                        int neg_or_pos = ((2 * (rand() % POSITIVE_FACTOR)) - 1);
                        element_value = random * neg_or_pos;
                    } else if (matrix_genr_state & INPUT_MODE) {
                        scanf("%lf", &element_value);
                    }
                    B[ i * N + j ] = element_value;
                }
            }
        }

        printf(" Finished populating Matricies.. This is the run..\n\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Finally perform the parallelization 
        of the matrix multiplciation */

    if (A != NULL && C != NULL) {

        double prev_time_sum = time_local;
        time_local = MPI_Wtime();

        // scatter into other processes
        err = MPI_Bcast(B, N*N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        // err = MPI_Bcast(A, N*N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        time_local = MPI_Wtime() - time_local;

        for (int i = 0; i < proc_count; ++i) {
            if (my_rank == i) {
                printf("My Rank: %d Broadcast took time=%lf\n", my_rank, time_local);
                printf("-----------------------\n");
                MPI_Barrier(MPI_COMM_WORLD);  
            } else {
                MPI_Barrier(MPI_COMM_WORLD); 
            }
        }

        prev_time_sum += time_local;

        time_local = MPI_Wtime();

        /* KIJ is a special case of scattering:
            need to scatter the vertical partitions of A...
            so that we have all K elements in A to perform a complete partial 
            dot product.. In order to scatter a vertical partition we need to 
            specify a stride value to MPI_Scatter in order to send partial 
            horizontal matricies of A */

            
        err = MPI_Scatter(A, N*rows_pproc, MPI_DOUBLE, A, N*rows_pproc, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

        // used in partial vector calculations
        register double sum = 0;
        //double part_sum = 0;
        double a_i_k = 0;  

        double* rowA;
        double* rowB;
        double* rowC;

        /* each of the cases is a matrix multiplication 
            implementation, of form ijk, ikj, kij */

        switch (matrix_mult_state){
            case NULL_MODE:
                printf("BAD MODE.. NOTHING DONE..\n");
            break;
            
            case IJK_MODE:
                /*  Modified to take into consideration cache miss rate on non contiguous memeory
                    when accessing multiple rows, this will likely result in many cache misses
                    I am now making sure that the next memory access is done on an address spacially 
                    close the to the previous address.. by iterating over a Transpose(B) */
                for (int i = 0; i < A_row_count; ++i) {
                    rowC = &(C[i*N]);
                    rowA = &(A[i * N]);
                    for (int j = 0; j < B_col_count; ++j) {
                        sum = 0;
                        rowB = &(B[j * N]);
                        for (int k = 0; k < vector_size; ++k) {
                            // slower by default since accessing distinct rows,
                            // even though data exists in linear contiguous array.. 
                            // sum += MatrixEntry(A,i,k,N) * MatrixEntry(B,k,j,N);

                            // faster access minimize cache 
                            // miss when accesing transposed matrix
                            // as fast as ikj form now..
                            sum += rowA[k] * rowB[k]; 
                        }
                        rowC[j] = sum;
                    }
                }
            break;
            
            case IKJ_MODE:
                /* for this mode, i optimized by reading in the second matrix 
                   in transposed form. Now iteration over the elements in a row
                   in B is faster. */
                for (int i = 0; i < A_row_count; ++i) {
                    rowC = &(C[i * N]);
                    rowA = &(A[i * N]); 
                    for (int k = 0; k < vector_size; ++k) {
                        a_i_k = rowA[k]; // MatrixEntry(A,i,k,N);
                        rowB = &(B[k*N]);
                        for (int j = 0; j < B_col_count; ++j) {
                            rowC[j] += a_i_k * rowB[j];
                        }
                    }
                }
            break;

            case KIJ_MODE:

                /*  considering cache efficiency by caching Aik and even the pointers 
                     to rows in B and C to access it, (hopefully) faster.

                    in the second loop i cache Aik but i think this could be sped up
                    even more by iterating over contiguous elements vs discondiguous elements

                    This can be improved by iterating over the elements in a row rather 
                    than over each row to access an element..

                    To do this improvement would require that the matrix A is sent transposed
                    when the matrix mult being requested is kij.. this implementation below
                    will need to be altered so that A is iterated over each element in a row
                    by iterating A[k*N + i].. but if and only if the A is transposed 

                    the parallel implementation will have to change... need to read in a 
                    transposed matrix, then scatter column groups (rather than row groups) */
                for (int k = 0; k < vector_size; ++k) {
                    rowB = &(B[k*N]);
                    rowA = &(A[k*A_row_count]);
                    for (int i = 0; i < A_row_count; ++i) {
                        rowC = &(C[i*N]);
                        a_i_k = rowA[i];    // transposed column major iteration (assumsed matrix is read in transposed..)
                        // a_i_k = A[i*N + k]; // standard matrix A_row_count tall
                        for (int j = 0; j < B_col_count; ++j) {
                            rowC[j] += a_i_k * rowB[j];
                        }
                    }
                }
            break;

            default:
            break;

        }

        rowA = NULL;
        rowB = NULL;
        rowC = NULL;

        if (my_rank == MASTER && rem_rows != 0) {
            printf("   Not Yet Handling case where N=%d rows are not evenly \n    distributable amongst p=%d processors.\n",  N, comm_sz);
            printf("   There are I=%d rows remaining.\n",  rem_rows);
        }

        err = MPI_Gather(  C, N*rows_pproc, MPI_DOUBLE,
                           C, N*rows_pproc, MPI_DOUBLE,
                           MASTER, MPI_COMM_WORLD);


        /* further implementation of calculating the residue matrix calculation.
            There are a remaining number of rows from A that have not yet been 
            evenly partitioned amongst the processes */

        if ( A_REM != NULL && 
             C_REM != NULL )
        {
            printf("Not Yet Implemented Yet! (TODO: calculate remaining rows)\n");
        }  


        time_local =  (MPI_Wtime() - time_local) + prev_time_sum;

        MPI_Reduce(&time_local, &elapsed, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    } else {
        /* this processor is not being utilized 
            since there are more PROCESSORS than ROWS in this case */
    }

    if (my_rank == MASTER) {
        
        printf("Matrix Multiplication Complete.\n elapsed time = %0.10f seconds\n",  elapsed);
        printf("\n");
        #if STDIN_RESULT_OUTPUT == true
        printf("Final Result:  \n");
        matrix_print(C, N, N);
        #endif
        printf("------------------------\n");
        printf("\n");
        
    }

    free(A);
    free(B);
    free(C);
    free(A_REM);
    free(C_REM);

    // err = MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
    //        MPI_Op op, int root, MPI_Comm comm)

    MPI_Finalize();

    return 0;
}


// for (int i = 0; i < proc_count; ++i) {
//     if (my_rank == i) {
//         printf("Proc Rank: #%d\n", my_rank);
//         printf("    Matrix: \n");
//         matrix_print(data->M, data->m, data->n);
//         MPI_Barrier(MPI_COMM_WORLD);  
//     } else {
//         MPI_Barrier(MPI_COMM_WORLD); 
//     }
// }
