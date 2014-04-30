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



#define MASTER 0

#define CHUNKSIZE   10

#define STDOUT_PRINT_RESULTS false

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
    TODO: parallelize this.. */
void gauss_reduce_parallel(double* A, int aug_row_size, int aug_col_size, int thread_count) {
    double* upper_triangle  = NULL;
    double* tmp_vector;

    /* paralellization values */
    const int chunk = CHUNKSIZE;
    int my_rank = 0;
    const int dh = aug_row_size/thread_count;
    int a = 0;
    int b = 0;

    double** B = (double**) malloc (sizeof(double*) * aug_row_size);

    // #pragma omp parallel
    for (int i = 0; i < aug_row_size; ++i) { 
        B[i] = &A[i*aug_row_size];
    }

    /* we wish to eliminate 'row' from the matrix and 
       make it's lower diagnal values equal 0 */
    for (int row = 1; row < aug_col_size; ++row) {
        for (int item = 0; item < row; ++item) {
            double  multiplier = B[row][item] 
                                    / B[item][item];
            
            #pragma omp parallel num_threads(thread_count) \
                        shared(B, multiplier) private(a,b,my_rank)
            my_rank = omp_get_thread_num();  //some unavoidable minimal overhead 
            a = my_rank * dh;
            b = a + dh;
            for (int i = a; i < b; ++i) {
                B[ row ][ i ] -= B[ item ][ i ] * multiplier;
            }
        }
    }

    free(B);
}

/* Perform gaussian elimination on the augmented matrix:
    A is an AUGMENTED MATRIX in this case..
    TODO: parallelize this.. */
double** gauss_reduce_partial_pivot_parallel(double* A,        int aug_row_size, 
                                             int aug_col_size, int thread_count) {
    double* upper_triangle  = NULL;
    double* tmp_vector;
    double  diag_element = 0;
    double  elim_element = 0;

    /* paralellization values */
    const int chunk = CHUNKSIZE;
    int my_rank = 0;
    const int dh = aug_row_size/thread_count;
    int a = 0;
    int b = 0;

    double** B = (double**) malloc (sizeof(double*) * aug_col_size);

        for (int i = 0; i < aug_col_size; ++i) {
            B[i] = &A[i*aug_row_size];
        }

        /* we wish to eliminate 'row' from the matrix and 
           make it's lower diagnal values equal 0 */

        /* Item is a column iteration that selects which column we are iterating over */
        /* Row will iterate over the rows in the item column.. */
        for (int item = 0; item < aug_row_size-1; ++item) {   
            for (int row = item+1; row < aug_col_size; ++row) { // row we with to eliminate
                /* for the rest of the values in the row up to the diagonal, 
                    need to access the mulipliers to eliminate the rest of the 
                    values in the row in the lower traingle (BELOW the diagonal) */

                // printf("[row][item]:[%d][%d] \n", row, item);

                elim_element = B[row][item];
                if (elim_element == 0) {
                    // then done.. no need for elimination its already eliminated
                    // printf("TODO: HANDLE CASE WHEN ROW TO ELIMINATE IS ALREADY ELIMINATED\n");
                } else {
                    /* First check for DIVIDE BY ZERO (which means the item to use is already eliminated
                        Therefore:, if 0 then either:
                            the item row should be in a new row in the matrix 
                            the item row should be in the same place, but if it is
                                then that measn the values before the 0 are also 0 since 
                                we are accessing that diag element because it is the first 
                                in the list..
                        Therefore: just swap */
                    diag_element = B[item][item];
                    if (diag_element == 0) {
                        // swap item row with 'row' row  
                        // (we wish to eliminate 'row' row)
                        // printf("TODO: HANDLE CASE WHEN DIAG HAPPENS TO BE ELIMINATED, SWAP WITH 'ROW'\n");
                        double* row_tmp_p = B[row];
                        B[row] = B[item];
                        B[item] = row_tmp_p;
                        row_tmp_p = NULL;
                    } else {

                        /* partial pivoting */
                        double max_value = fabs(B[row][item]);
                        // printf("[row][item]:[%d][%d]Max Value = %lf\n", row, item, max_value);
                        int    max_val_row = row;
                        for (int i = a+1; i < aug_col_size; ++i) {
                            if (max_value < fabs(B[i][item])) {
                                max_value = fabs(B[i][item]);
                                max_val_row = i;
                            }
                        }

                        if (max_value > fabs(B[row][item])) {
                            // swap the max_val_row row with the 'row' row
                            // matrix_row_swap(A)             // to much overhead..
                            double* row_tmp_p = B[max_val_row]; // just swap pointers..
                            B[max_val_row]    = B[row];
                            B[row] = row_tmp_p;
                            row_tmp_p = NULL;

                            elim_element = B[row][item];
                            diag_element = B[item][item];
                        } // else no swap necessary
                    

                        /* perform gaussian elimination on the row
                        Row subraaction to reduce a row */
                        #pragma omp parallel num_threads(thread_count) \
                                    shared(B) private(a,b,my_rank)
                        {
                            my_rank = omp_get_thread_num();  //some unavoidable minimal overhead 
                                                    /* find the multiplier factor */
                            double  multiplier = elim_element / diag_element;
                            a = my_rank * dh;
                            b = a + dh;
                            // #pragma omp for 
                            for (int i = a; i < b; ++i) {
                                B[ row ][ i ] -= B[ item ][ i ] * multiplier;
                            }
                        }
                        // end pragma parallel region


                    } //end if 0'd diag
                } // end element already eliminated

            } // end for each row from 1+diagonal to end of column size
        } // end for each column except for the agumented column

    
    return B;
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
            /* perform substitution on each col item wit the current x
                then reduce it a bit so that it zero's out the matrix's 
                element since we subtract both sides by some partialy evaluated 
                term at VALUE x   */
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

#if STDOUT_PRINT_RESULTS == true
    matrix_print(augmented_matrix, N+1, N);
#endif
    printf("\n");

    double** reduced_echelon_2d; 
    double*  reduced_echelon_1d;

    /* Serial Gaussian Elimination  no pivot */

    double* B2      = matrix_copy(augmented_matrix, aug_row_size * aug_col_size);
        /* find the upper triangle */
    start_time      = omp_get_wtime();
    gauss_reduce_serial(B2, aug_row_size, aug_col_size);
    start_time      = omp_get_wtime() - start_time;
    printf("Gaussian Serial =                                       %lf\n", start_time);

#if STDOUT_PRINT_RESULTS == true
    matrix_print(B2, aug_row_size, aug_col_size);
#endif
/* BACK SUBSTITION ON 1D MATRIX */
    start_time                  = omp_get_wtime();
    reduced_echelon_1d          = row_reduce(B2, aug_row_size, aug_col_size);
    start_time                  = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_RESULTS == true
    matrix_print(reduced_echelon_1d, N+1, N);
#endif
    printf("x       =: ");
    x_result = matrix_col_get(reduced_echelon_1d, N, N+1, N);
    vector_print( x_result, N );


    row_count = N;  // row size 
    col_count = N; // column size
    b_result  = matrix_mult_vector(A, x_result, row_count, col_count);
    printf("i2norm  =: %lf\n", LSquareNorm( vector_subtract(b_result, b, row_count, NULL), row_count) );

    printf("\n");


    /* Parallel Gaussian Elimination  */

    double* copy_aug_mat2   = matrix_copy(augmented_matrix, aug_row_size * aug_col_size);

    /* find the upper triange */
    start_time      = omp_get_wtime();
    gauss_reduce_parallel(copy_aug_mat2, aug_row_size, aug_col_size, number_of_threads);
    start_time      = omp_get_wtime() - start_time;
    printf("Gaussian Parallel Time=                                 %lf\n", start_time);

#if STDOUT_PRINT_RESULTS == true
    matrix_print(copy_aug_mat2, N+1, N);
#endif

    /* BACK SUBSTITION ON 1D MATRIX */
    start_time                  = omp_get_wtime();
    reduced_echelon_1d          = row_reduce(copy_aug_mat2, aug_row_size, aug_col_size);
    start_time                  = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_RESULTS == true
    matrix_print(reduced_echelon_1d, N+1, N);
#endif
    printf("x       =: ");
    x_result = matrix_col_get(reduced_echelon_1d, N, N+1, N);
    vector_print( x_result, N );

    row_count = N;  // row size 
    col_count = N; // column size
    b_result  = matrix_mult_vector(A, x_result, row_count, col_count);
    printf("i2norm  =: %lf\n", LSquareNorm( vector_subtract(b_result, b, row_count, NULL), row_count) );

printf("\n");


    /* Gaussian Elimination w/ Partial Pivoting Test (PARALLEL) */

    double* copy_aug_mat   = matrix_copy(augmented_matrix, aug_row_size * aug_col_size);

    /* find the upper triangle using partial pivoting */
    start_time = omp_get_wtime();
    double** B = gauss_reduce_partial_pivot_parallel(copy_aug_mat, aug_row_size, aug_col_size, number_of_threads);
    start_time = omp_get_wtime() - start_time;
    printf("Gaussian Parallel w/ Partial Pivot. Time =              %lf\n", start_time);

#if STDOUT_PRINT_RESULTS == true
    // matrix_print(copy_aug_mat, N+1, N);
    matrix_2d_print(B, aug_row_size, aug_col_size);
#endif

    /* BACK SUBSTITION ON 2D MATRIX */
    start_time                  = omp_get_wtime();
    reduced_echelon_2d          = row_reduce_2d(B, aug_row_size, aug_col_size);
    start_time                  = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_RESULTS == true
    matrix_2d_print(reduced_echelon_2d, N+1, N);
#endif
    printf("x       =: ");
    x_result = matrix_col_get(reduced_echelon_2d[0], N, N+1, N);
    vector_print( x_result, N );


    row_count = N;  // row size 
    col_count = N; // column size
    b_result  = matrix_mult_vector(A, x_result, row_count, col_count);
    printf("i2norm  =: %lf\n", LSquareNorm( vector_subtract(b_result, b, row_count, NULL), row_count) );

printf("\n");
printf("\n");

    matrix_delete(A);
    matrix_delete(reduced_echelon_1d);
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