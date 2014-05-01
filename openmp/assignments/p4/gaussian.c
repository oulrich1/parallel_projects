 /* 
    Gaussian Elimincation with Partial Pivoting
        in Parallel:

    By Oriah Ulrich

    Solves Ax = b
    by performing gaussian elimination on A in order to find x

    outputs the time it takes to perform some of the major operations
        - generate matrix
        - gaussian serial
        - gaussian parallel
        - gassian parallel with rowwise partial pivoting
        - row reduced echelon 

    outputs the result of the operations:
        - x resulting vector
        - i2norm of the residual differences between 
            b_result and given b  (ie: Ax - b)
    
    host:juch-s04@207.108.8.131 
    pass:!Q2..
*/

    /* MUST COMPILE WITH -lm */



#include "matrix.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
    // #include <pthread.h>

#ifndef NULL
#define NULL 0
#endif

typedef unsigned int uint;
typedef long double LD_T;

typedef double PRIM_T;


#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define equ(a, b) ((a) == (b) ? (true) : (false))

#define MASTER 0

#define CHUNKSIZE   10


#ifndef TEST_SERIAL_GAUSS
#define TEST_SERIAL_GAUSS false
#endif

#ifndef TEST_PARALLEL_GAUSS
#define TEST_PARALLEL_GAUSS true
#endif

#ifndef STDOUT_PRINT_MATRIX_RESULTS
#define STDOUT_PRINT_MATRIX_RESULTS false
#endif

#ifndef STDOUT_PRINT_X_RESULTS
#define STDOUT_PRINT_X_RESULTS false
#endif

/* Perform gaussian elimination on the augmented matrix:
    A is an AUGMENTED MATRIX in this case..
    TODO: parallelize this.. */
double* gauss_reduce_serial(double* A, int aug_row_size, int aug_col_size) {
    double* upper_triangle  = NULL;
    double* tmp_vector;

    /* we wish to eliminate 'row' from the matrix and 
       make it's lower diagnal values equal 0 */
    for (int row = 1; row < aug_col_size; ++row) {
        int elements_skip_eliminate_iter = row * aug_row_size;
        int elements_skip_factor_iter = 0;
        for (int item = 0; item < row; ++item) {

            elements_skip_factor_iter = item * aug_row_size;
            double  multiplier = A[elements_skip_eliminate_iter + item] 
                                    / A[elements_skip_factor_iter + item];

            for (int i = 0; i < aug_row_size; ++i) {
                A[ elements_skip_eliminate_iter + i ] -= A[ elements_skip_factor_iter + i ] * multiplier;
            }
        }
    }

    return A;
}

/* Perform gaussian elimination on the augmented matrix:
    A is an AUGMENTED MATRIX in this case..
    TODO: parallelize this.. returns A */
double** gauss_reduce_partial_pivot_serial(double** A,       int aug_row_size, 
                                             int aug_col_size, int thread_count) {
    double* upper_triangle  = NULL;
    double* tmp_vector;
    double  diag_element = 0; // is the greatest value by the end of partial pivoting step
    double  elim_element = 0; // is the element in the column vector below diagonal that 
                              //    is to be eliminated
    double  multiplier = 1;

    double  greatest_value = 0;
    int     greatest_val_row_id = 0;

    /* iterators */
    int row = 0,
        pivot = 0;
    int k = 0;

    /* iterate through the diagonal */
    for (pivot = 0; pivot < min(aug_col_size, aug_row_size); ++pivot) { 

        /*  FIND GREATEST VALUE IN COLUMN from the 
            pivot position element in the row..  */
        greatest_value = A[pivot][pivot];
        greatest_val_row_id = pivot;

        for (row = pivot + 1; row < aug_col_size; ++row) {
            if (labs(greatest_value) < labs(A[row][pivot])) { 
                greatest_value = A[row][pivot];
                greatest_val_row_id = row;
            }
        }

        if(greatest_value == 0) {
            printf("ERROR: Pivot was found to be 0.. this system is not solvable with only gaussian elimination \n");
            return NULL;
        }

        /* swap */
        if (greatest_val_row_id != pivot) { // then swap
            double* tmp_greatest_coef_row_ptr = A[greatest_val_row_id];
            A[greatest_val_row_id] = A[pivot];
            A[pivot] = tmp_greatest_coef_row_ptr;
            tmp_greatest_coef_row_ptr = NULL;
        }

        /* for the rest of the rows down from the pivot position row
            execute a gaussian elimination step for each row using the 
            pivoted row as the muliplier */

        for (row = pivot + 1; row < aug_col_size; ++row) {
            multiplier = A[row][pivot] / A[pivot][pivot];

            for (k = 0; k < aug_row_size; ++k) {
                A[ row ][ k ] -= A[ pivot ][ k ] * multiplier;
            } 
        }  
    } // end for all elements in diagonal

    return A;
}

/* Perform gaussian elimination on the augmented matrix:
    A is an AUGMENTED MATRIX in this case..
    TODO: parallelize this.. returns A */
double** gauss_reduce_partial_pivot_parallel(double** A,       int aug_row_size, 
                                             int aug_col_size, int thread_count) {
    double* upper_triangle  = NULL;
    double* tmp_vector;
    double  diag_element = 0; // is the greatest value by the end of partial pivoting step
    double  elim_element = 0; // is the element in the column vector below diagonal that 
                              //    is to be eliminated

    /* paralellization values */
    const int chunk = CHUNKSIZE;
    double  multiplier = 1;
    int     my_rank = 0;

    const int dh = aug_row_size/thread_count;
    int     a = 0;
    int     b = 0;

    double  greatest_value = 0;
    int     greatest_val_row_id = 0;

    double  local_great_val = 0;
    int     local_great_val_row_id = 0;

    /* iterators */
    int row = 0,
        pivot = 0;
    int k = 0;

    /* iterate through the diagonal */

    for (pivot = 0; pivot < min(aug_col_size, aug_row_size); ++pivot) { 

        /*  FIND GREATEST VALUE IN COLUMN from the 
            pivot position element in the row..  */
        greatest_value = A[pivot][pivot];
        greatest_val_row_id = pivot;

        for (row = pivot + 1; row < aug_col_size; ++row) {
            if (labs(greatest_value) < labs(A[row][pivot])) { 
                greatest_value = A[row][pivot];
                greatest_val_row_id = row;
            }
        }

        if(greatest_value == 0) {
            printf("ERROR: Pivot was found to be 0.. this system is not solvable with only gaussian elimination \n");
            // return NULL;
        }

        /* swap */
        if (greatest_val_row_id != pivot) { // then swap
            double* tmp_greatest_coef_row_ptr = A[greatest_val_row_id];
            A[greatest_val_row_id] = A[pivot];
            A[pivot] = tmp_greatest_coef_row_ptr;
            tmp_greatest_coef_row_ptr = NULL;
        }

        /* for the rest of the rows down from the pivot position row
            execute a gaussian elimination step for each row using the 
            pivoted row as the muliplier */

        #pragma omp parallel num_threads(thread_count) \
            shared(A, pivot, aug_col_size) \
            private(row, multiplier, k)
        {
            #pragma omp for
            for (row = pivot + 1; row < aug_col_size; ++row) {
                multiplier = A[row][pivot] / A[pivot][pivot];
                /* perform a standard row elimination step */
                for (k = 0; k < aug_row_size; ++k) {
                    A[ row ][ k ] -= A[ pivot ][ k ] * multiplier;
                } 
            }

        } // end pragma

        
    }

    return A;
}

            //vector_scalar_mult( , aug_col_size, x, NULL)

/* performs a backward substitution row reduction */
double* row_reduce(double* A, int aug_row_size, int aug_col_size){
    double* reduced_echelon = NULL;
    int b_vec_id   = aug_row_size - 1;
    int A_row_size = aug_row_size - 2;

    double  x       = 0;
    int     b_iter  = 0;
    int     diag_iter = 0;
    int     col_iter = 0; 
    for (int row = aug_col_size - 1; row >= 0; --row)
    {
        b_iter = row * aug_row_size + b_vec_id;
        diag_iter = row* aug_row_size + row;
        x = A[ b_iter ] / A[ diag_iter ];
        A[ b_iter ] = x;
        A[ diag_iter ] = 1;

        col_iter = 0;
        for (int col_item = row - 1; col_item >= 0; --col_item)
        {
            b_iter   = col_item * aug_row_size + b_vec_id;
            col_iter = col_item * aug_row_size + row;
            A[ b_iter ] -= ( A[ col_iter ] * x );
            A[ col_iter ] = 0;
        }
    }

    return A; 
}


/* performs a backward substitution row reduction */
double** row_reduce_2d(double** B, int aug_row_size, int aug_col_size){
    double* reduced_echelon = NULL;
    int b_vec_id   = aug_row_size - 1;
    int A_row_size = aug_row_size - 2;
    double x = 0;
    for (int row = aug_col_size - 1; row >= 0; --row) {
        x = B[row][b_vec_id] / B[row][row];
        B[row][b_vec_id] = x;
        B[row][row] = 1;
        for (int col_item = row - 1; col_item >= 0; --col_item) {
            B[col_item][b_vec_id] -= ( B[col_item][row] * x );
            B[col_item][row] = 0;
        }
    }
    return B; 
}

double LSquareNorm(double* v1, int size) {
    double sum_squares = 0;
    for (int i = 0; i < size; ++i) {
        sum_squares += v1[i] * v1[i];
    }
    return sqrt(sum_squares);
}

void usage(){
    printf("\n");
    printf(" Usage: ./gaussian n p\n");
    printf("\n");
    printf("        n - size of the matrix [n x n]\n");
    printf("        p - number of threads to utilize\n");
    printf("\n");
}


/*  -----------------------------------------------------------------

    1.) Read in arguments
    2.) Randomly generate a matrix A and a column vector b
    3.) Gaussian elimination on Ax = b to find x
    4.) Using the resulting x, double check we got the right x_result
        by calculating b_result and finding the l2norm of the b's. 

*/
int main(int argc, char const *argv[])
{       
    if (argc != 3) {
        usage();
        exit(1);
    }

    srand(time(NULL));

    int N = atoi(argv[1]);  // size if matrix and b vector
    int P = atoi(argv[2]);  // number of threads for openmp

    int number_of_cores = omp_get_num_procs();
    int number_of_threads =  P; // omp_get_num_threads();

    printf("number of cores   = %d\n", number_of_cores);
    printf("number of threads = %d\n", number_of_threads);

    int A_row_size = N;
    int A_col_size = N;
    int b_size = N;

    /* not available */
    // struct rusage usage;
    // getrusage(P, &rusage);

    double* A = matrix_new(N, N);
    double* b = vector_new(N);

    /* used for calculating the i2norm of the residual */
    int     row_count = N,  // row size 
            col_count = N;  // column size

    double* x_result = NULL; // x result from our gaussian elimination
    double* b_result = NULL; // A * x_result

    int aug_row_size = N+1;
    int aug_col_size = N;
    double* augmented_matrix = matrix_new(aug_col_size, aug_row_size);

    double start_time = omp_get_wtime();
    vector_randomize(b, b_size);
    matrix_randomize(A, A_row_size, A_col_size);
    matrix_copy_columns (A, A_row_size, A_col_size,
                         augmented_matrix, aug_row_size);
    /* augment matrix A with a new column vector b .. 
        for a consistent gaussian experience.. */
    matrix_col_set( augmented_matrix, aug_row_size-1, aug_row_size, 
                    b, b_size);
    start_time = omp_get_wtime() - start_time;
    printf("Generated Matrix:                                       %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_print(augmented_matrix, N+1, N);
#endif
    printf("----------------------------\n");


    double** reduced_echelon_2d; 
    double*  reduced_echelon_1d;

#if TEST_SERIAL_GAUSS == true
    /* Serial Gaussian Elimination  no pivot */

    double** B2     = matrix_copy_1d2d(augmented_matrix, aug_row_size,  aug_col_size);
        /* find the upper triangle */
    start_time      = omp_get_wtime();
    gauss_reduce_partial_pivot_serial(B2, aug_row_size, aug_col_size, 1);
    start_time      = omp_get_wtime() - start_time;
    printf("Gaussian Serial =                                       %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_2d_print(B2, aug_row_size, aug_col_size);
#endif
/* BACK SUBSTITION ON 1D MATRIX */
    start_time                  = omp_get_wtime();
    reduced_echelon_2d          = row_reduce_2d(B2, aug_row_size, aug_col_size);
    start_time                  = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_2d_print(reduced_echelon_2d, N+1, N);
#endif

    x_result = matrix2d_col_get(reduced_echelon_2d, N, N+1, N);
#if STDOUT_PRINT_X_RESULTS == true
    printf("x       =: ");
    vector_print( x_result, N );
#endif


    row_count = N;  // row size 
    col_count = N; // column size
    b_result  = matrix_mult_vector(A, x_result, row_count, col_count);
    printf("i2norm  =: %lf\n", LSquareNorm( vector_subtract(b_result, b, row_count, NULL), row_count) );

    printf("----------------------------\n");
#endif



#if TEST_SERIAL_GAUSS == true

    /* Serial Gaussian Elimination, without partial pivoting  */

    double* copy_aug_mat2   = matrix_copy(augmented_matrix, aug_row_size * aug_col_size);

    /* find the upper triange */
    start_time      = omp_get_wtime();
    gauss_reduce_serial(copy_aug_mat2, aug_row_size, aug_col_size);
    start_time      = omp_get_wtime() - start_time;
    printf("Gaussian Parallel Time=                                 %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_print(copy_aug_mat2, N+1, N);
#endif

    /* BACK SUBSTITION ON 1D MATRIX */
    start_time                  = omp_get_wtime();
    reduced_echelon_1d          = row_reduce(copy_aug_mat2, aug_row_size, aug_col_size);
    start_time                  = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_print(reduced_echelon_1d, N+1, N);
#endif

    x_result = matrix_col_get(reduced_echelon_1d, N, N+1, N);
#if STDOUT_PRINT_X_RESULTS == true
    printf("x       =: ");
    vector_print( x_result, N );
#endif

    row_count = N;  // row size 
    col_count = N; // column size
    b_result  = matrix_mult_vector(A, x_result, row_count, col_count);
    printf("i2norm  =: %lf\n", LSquareNorm( vector_subtract(b_result, b, row_count, NULL), row_count) );

printf("----------------------------\n");

#endif

    /* Gaussian Elimination w/ Partial Pivoting Test (PARALLEL) */

#if TEST_PARALLEL_GAUSS == true

    double** copy_aug_mat   = matrix_copy_1d2d(augmented_matrix, aug_row_size, aug_col_size);

    /* find the upper triangle using partial pivoting */
    start_time = omp_get_wtime();
    double** B = gauss_reduce_partial_pivot_parallel(copy_aug_mat, aug_row_size, aug_col_size, number_of_threads);
    start_time = omp_get_wtime() - start_time;
    printf("Gaussian Parallel w/ Partial Pivot. Time =              %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    // matrix_print(copy_aug_mat, N+1, N);
    matrix_2d_print(B, aug_row_size, aug_col_size);
#endif

    /* BACK SUBSTITION ON 2D MATRIX */
    start_time                  = omp_get_wtime();
    reduced_echelon_2d          = row_reduce_2d(B, aug_row_size, aug_col_size);
    start_time                  = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_2d_print(reduced_echelon_2d, N+1, N);
#endif

    x_result = matrix2d_col_get(reduced_echelon_2d, N, N+1, N);
#if STDOUT_PRINT_X_RESULTS == true
    printf("x       =: ");
    vector_print( x_result, N );
#endif

    row_count = N;  // row size 
    col_count = N; // column size
    b_result  = matrix_mult_vector(A, x_result, row_count, col_count);
    printf("i2norm  =: %lf\n", LSquareNorm( vector_subtract(b_result, b, row_count, NULL), row_count) );

printf("----------------------------\n");
printf("\n");

#endif

    matrix_delete(A);
    // matrix_delete(reduced_echelon_1d);
    matrix_2d_delete(reduced_echelon_2d);
    matrix_delete(augmented_matrix);
    vector_delete(b);

    return 0;
}


// #pragma omp parallel for num_threads(thread_count)
                        // parallel for shared(B) 


// #pragma omp parallel num_threads(thread_count) private(my_rank)

// #pragma omp for schedule(dynamic, 1000) nowait
// for (int i = 0; i < aug_row_size; ++i) {
//     B[ item ][ i ] *= multiplier;
// }

// parallel for shared(B)
// #pragma omp section 
// #pragma omp for schedule(dynamic, 1000) nowait


        /* tried to reduce on MAX value */
        // #pragma omp parallel shared(A, row, pivot, greatest_value, greatest_val_row_id) private(local_great_val, local_great_val_row_id)
        // {
        //     #pragma omp for
        //     for (row = pivot + 1; row < aug_col_size; ++row) {
        //         if (labs(local_great_val) < labs(A[row][pivot])) { 
        //             local_great_val = A[row][pivot];
        //             local_great_val_row_id = row;
        //         }
        //     }

        //     #pragma omp critical 
        //     {   if (labs(local_great_val) > labs(greatest_value)) {   
        //             greatest_value = local_great_val;
        //             greatest_val_row_id = local_great_val_row_id;
        //         }
        //     }
        // }
