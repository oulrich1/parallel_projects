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



// #include "matrix.h"

#include <time.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
// #include <unistd.h> /* machine info */

#include <string.h>
#include <math.h>

#include <omp.h>

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

#ifndef MAX_VECTOR_VALUE
#define MAX_VECTOR_VALUE 1000
#endif

double drand() {
    return (rand()/(double)RAND_MAX); // c99 doesnt support drand, generate 0 to 1 value
}

double randomDouble() {
    return ((2 * (rand() % 2) - 1) * drand() * (MAX_VECTOR_VALUE));
}

/* Matrix Definitions */

typedef struct Matrix
{
    double** data;
    int N;
    int M;
#define Entry(MAT,i,j) (MAT->data[i][j])
} Matrix;

/* Just a Matrix but with helpers to access specifics */
typedef struct MatrixAugmented
{
    Matrix mat;
    int b_vec_index;
#define MAT(mm) ((mm)->mat)
#define MATDATA(mm) (((mm)->mat).data)
} MatrixAugmented;

/* ---------------- Matrix Method Definitions ---------------- */

Matrix* matrix_new(){
    Matrix* matrix = (Matrix*) malloc(sizeof(Matrix));
    matrix->M = 0;
    matrix->N = 0;
    matrix->data = NULL; // NULL
    return matrix;
}

/* must be freed */
void matrix_free(Matrix* m1){
    if (!m1) {
        return;
    }
    for (int i = 0; i < m1->M; ++i) {
        free(m1->data[i]);
    }
    free(m1->data);
}

void matrix_init(Matrix* mat, int M, int N) {
    if (!mat) {
        return;
    }
    mat->M = M;
    mat->N = N;
    mat->data = (double**) malloc (sizeof(double*) * M);
    for (int i = 0; i < M; ++i) {
        mat->data[i] = (double*) calloc (N, sizeof(double));
    }
}

void matrix_randomize(Matrix* mat) {
    if (!mat || !mat->data ) {
        return;
    }
    for (int i = 0; i < mat->M; ++i) {
        for (int j = 0; j < mat->N; ++j) {
            mat->data[i][j] = randomDouble();
        }
    }
}

void matrix_resize(Matrix* mat, int M, int N) {
    if (!mat) {
        return;
    }
    matrix_init(mat, M, N);
}

void matrix_print(Matrix* mat) {
    if (!mat) {
        return;
    }
    int col_size = mat->M; // # rows
    int row_size = mat->N; // # cols
    printf("[ \n");
    for (int i = 0; i < col_size; ++i)
    {
        printf("    [ ");
        for (int j = 0; j <= row_size-2; ++j)
        {
            printf("%.2lf, ", Entry(mat,i,j) );
        }
        printf("%.2lf ]", Entry(mat,i,row_size-1) );
        if (i != col_size-1) {
            printf(",\n");
        } else {
            printf("\n");
        }
    }
    printf("] \n\n");
}

/* Augment Matrix Method to augment two matricies m1 and m2 into aug */

void augmat_augment(MatrixAugmented *aug, Matrix *m1, Matrix *m2) {
    if (!aug || !m1 || !m2) {
        return;
    }
    aug->b_vec_index = m1->N;
    matrix_resize(&(aug->mat), m1->M, m1->N + m2->N);
    for (int i = 0; i < m1->M; ++i) {
       for (int j = 0; j < m1->N; ++j) {
            aug->mat.data[i][j] = m1->data[i][j];         
       }
    }
    int aug_N = aug->mat.N-1;
    for (int i = 0; i < m2->M; ++i) {
        for (int j = 0; j < m2->N; ++j) {
            aug->mat.data[i][aug_N+j] = m2->data[i][j];
        }
    }

}

Matrix* matrix_col_get (Matrix* m, int index) {
    if (!m) {
        return NULL;
    }
    Matrix* matrix = matrix_new();
    matrix_init(matrix, 1, m->M); // 1 x M row vector.. for efficiency
    for (int i = 0; i < m->M; ++i) {
        matrix->data[0][i] = m->data[i][index];
    }
    return matrix;
}


/* assumes col count is the size of matrix,
    and matrix is square and that x_result must be a vector
    returns a new matrix resulting from a matrix multiply */
Matrix* matrix_mult(Matrix* m1, Matrix* m2) {
    Matrix* C = matrix_new();
    matrix_init(C, m1->M, m2->N);
    if (m1->N != m2->M) {
        printf("m1 row size not equal to m2 col size..\n");
    }

    double sum = 0;
    for (int i = 0; i < m1->M; ++i) {
        for (int j = 0; j < m2->N; ++j) {
            sum = 0;
            for (int k = 0; k < C->M; ++k) {
                sum += m1->data[i][k] * m2->data[k][j];
            }
            C->data[i][j] = sum;
        }
    }

    return C;
}

/*  Subtracts a matrix m2 from m1: 
    returns a new matrix result */

Matrix* matrix_subtract(Matrix* m1, Matrix* m2) {
    Matrix* m = matrix_new();
    matrix_init(m, m1->M, m1->N);
    for (int i = 0; i < m1->M; ++i) {
        for (int j = 0; j < m1->N; ++j) {
            m->data[i][j] = m1->data[i][j] - m2->data[i][j];
        }
    }
    return m;
}

/*  Transposes a matrix: 
    returns a new matrix */

Matrix* matrix_transpose(Matrix* m1){
    if (!m1) {
        return NULL;
    }
    int M;
    int N;
    M = m1->N;
    N = m1->M;
    Matrix* transposed = matrix_new();
    matrix_init(transposed, M, N);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            transposed->data[i][j] = m1->data[j][i];
        }
    }
    return transposed;
}

/* ---------------- End matrix helpers ---------------- */


/* Perform gaussian elimination with partial 
    pivoting on the augmented matrix */
void gauss_reduce_partial_pivot_parallel(MatrixAugmented *augmat, int thread_count) {

    Matrix* mat = &(augmat->mat);

    int aug_col_size = mat->M;
    int aug_row_size = mat->N;

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
        greatest_value = mat->data[pivot][pivot];
    
        greatest_val_row_id = pivot;

        for (row = pivot + 1; row < aug_col_size; ++row) {
            if (labs(greatest_value) < labs(mat->data[row][pivot])) { 
                greatest_value = mat->data[row][pivot];
                greatest_val_row_id = row;
            }
        }
        if(greatest_value == 0) {
            printf("ERROR: Pivot was found to be 0.. this system is not solvable with only gaussian elimination \n");
            // return NULL;
        }
        /* swap */
        if (greatest_val_row_id != pivot) { // then swap
            double* tmp_greatest_coef_row_ptr = mat->data[greatest_val_row_id];
            mat->data[greatest_val_row_id] = mat->data[pivot];
            mat->data[pivot] = tmp_greatest_coef_row_ptr;
            tmp_greatest_coef_row_ptr = NULL;
        }
        /* for the rest of the rows down from the pivot position row
            execute a gaussian elimination step for each row using the 
            pivoted row as the muliplier */
        #pragma omp parallel num_threads(thread_count) \
            shared(mat, pivot, aug_col_size) \
            private(row, multiplier, k)
        {
        
            #pragma omp for
            for (row = pivot + 1; row < aug_col_size; ++row) {
                multiplier = mat->data[row][pivot] / mat->data[pivot][pivot];
                /* perform a standard row elimination step */
                for (k = 0; k < aug_row_size; ++k) {
                    mat->data[ row ][ k ] -= mat->data[ pivot ][ k ] * multiplier;
                } 
            }

        } // end pragma
        
    }
}

/* performs a backward substitution row reduction
   substitutes the entire column with the value that 
   has been 'discovered' */
void row_reduce(MatrixAugmented* augmat){
    // double* reduced_echelon = NULL;
    int b_vec_index = augmat->b_vec_index;
    // int A_row_size = MAT(augmat).N - 2;
    int aug_col_size = MAT(augmat).M;
    double x = 0;
    for (int row = aug_col_size - 1; row >= 0; --row) {
        x = MATDATA(augmat)[row][b_vec_index] / MATDATA(augmat)[row][row];
        MATDATA(augmat)[row][b_vec_index] = x;
        MATDATA(augmat)[row][row] = 1;
        for (int col_item = row - 1; col_item >= 0; --col_item) {
            MATDATA(augmat)[col_item][b_vec_index] -= ( MATDATA(augmat)[col_item][row] * x );
            MATDATA(augmat)[col_item][row] = 0;
        }
    }
}

/* magnitude of difference / error */
double LSquareNorm(Matrix* m1) {
    double sum_squares = 0;
    for (int i = 0; i < m1->M; ++i) {
        for (int j = 0; j < m1->N; ++j) {
            sum_squares += m1->data[i][j] * m1->data[i][j];
        }
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


    Matrix mat; // Ax = b  // this is A
    Matrix b;
    matrix_init(&mat, N, N);
    matrix_init(&b, N, 1);
    matrix_randomize(&mat);
    matrix_randomize(&b);

    MatrixAugmented augmat;
    augmat_augment(&augmat, &mat, &b);

#if TEST_PARALLEL_GAUSS == true

    /* find the upper triangle using partial pivoting */
    double start_time = omp_get_wtime();
    gauss_reduce_partial_pivot_parallel(&augmat, number_of_threads);
    start_time = omp_get_wtime() - start_time;
    printf("Gaussian Parallel w/ Partial Pivot. Time =   %lf\n", start_time);


#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_print(&augmat.mat);
#endif

    /* BACK SUBSTITION ON 2D MATRIX */
    start_time = omp_get_wtime();
    row_reduce(&augmat);
    start_time = omp_get_wtime() - start_time;
    printf("Row Reduced Echelon (Serial) Time:                      %lf\n", start_time);

#if STDOUT_PRINT_MATRIX_RESULTS == true
    matrix_print(&augmat.mat);
#endif

 // matrix_print(&augmat.mat);
    Matrix* x_result = matrix_col_get(&augmat.mat, augmat.mat.N-1);

#if STDOUT_PRINT_X_RESULTS == true
    printf("x       =: ");
    matrix_print( x_result );
#endif

    Matrix* x_transposed = matrix_transpose(x_result);
    Matrix* b_result  = matrix_mult(&mat, x_transposed);
    Matrix* differences = matrix_subtract(b_result, &b);
    printf("i2norm  =: %lf\n", LSquareNorm( differences));

    free(b_result);
    free(x_transposed);
    free(x_result);
    free(differences);
    matrix_free(&b);
    matrix_free(&mat);
    matrix_free(&augmat.mat);
printf("----------------------------\n");
printf("\n");

#endif

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
