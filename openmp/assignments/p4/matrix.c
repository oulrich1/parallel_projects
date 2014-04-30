#include "matrix.h"

inline double drand() {
    return (rand()/(double)RAND_MAX); // c99 doesnt support drand, generate 0 to 1 value
}

inline double randomDouble() {
    return ((2 * (rand() % 2) - 1) * drand() * (MAX_VECTOR_VALUE));
}

/* Some basic matrix method definitons */

double* matrix_new(int m, int n) {
    if (m <= 0 || n <= 0) {
        return NULL;
    } 
    return (double*) calloc(m*n, sizeof(double));
}

void matrix_randomize (double* A, int row_size, int col_size){
    int num_elements = row_size * col_size;
    for (int i = 0; i < num_elements; ++i) {
         A[i] = randomDouble();
    }
}

bool matrix_delete(double* matrix) {
    if (matrix) {
        free(matrix);
    }
    return true;
}

bool matrix_2d_delete(double** matrix2d){
    if(matrix2d){
        free(matrix2d[0]);
    }
    return true;
}

/* returns a pointer to a full copy of A */
double* matrix_copy(double* A, int size){
    double* B = (double*) calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i) {
         B[i] = A[i]; // [m*n] => size values set
    }
    return B;
}

/* A size = row_size * col_size must be smaller than or equal to 
   the matrix we wish to copy the data to.. column major order */
void matrix_copy_columns(double* A, int row_size,   int col_size, 
                         double* B, int b_row_size) 
{
    for (int j = 0; j < row_size; ++j) {
        for (int i = 0; i < col_size; ++i) {
            B[ i*b_row_size + j ] = A[ i*row_size + j ];
        }
    }
}

/* Performs a Set Data operation on a 
   contiguious region of memory */

/* assumes row major matrix, sets the data in 
   the row_id-th vector in A to row_vector */
bool matrix_row_set(double* A,          int row_id, 
                    double* row_vector, int row_size) {
    if (row_size <= 0) {
        return false;
    }
    for (int j = 0; j < row_size; ++j) {
        A[ row_id * row_size + j ] = row_vector[j];
    }
    return true;
}


/* assumes row major matrix, sets the data in 
   the col_id-th column vector in A to col_vector */
bool matrix_col_set(double* A,          int col_id,  int row_size,
                    double* col_vector, int col_size) {
    if (col_size <= 0 ||    row_size <= 0) {
        return false;
    }
    for (int i = 0; i < col_size; ++i) {
        A[i * row_size + col_id] = col_vector[i];
    }
    return true;
}

/*  Performs a Get Data operation from a contiguious region of memory 

    The returned vector must be deleted manually..
*/

/* assumes row major matrix, sets the data in 
   the row_id-th vector in A to row_vector */
double* matrix_row_get(double* A, int row_id, int row_size) {
    // if (row_id < 0 || row_size <= 0) {
    //     return NULL;
    // }
    double* row_vector = vector_new(row_size);
    for (int j = 0; j < row_size; ++j) {
        row_vector[j] = A[ row_id * row_size + j ];
    }
    return row_vector;
}

/* assumes row major matrix,  */
double* matrix_col_get(double* A, int col_id, int row_size, int col_size) {
    if (col_id < 0 || row_size <= 0 || col_size <= 0){
        return NULL;
    }
    double* col_vector = vector_new(col_size);
    for (int i = 0; i < col_size; ++i) {
        col_vector[i] = A[i * row_size + col_id];
    }
    return col_vector;
}

/* assumes col count is the size of matrix,
    and matrix is square and that x_result must be a vector
    returns a new matrix resulting from a matrix multiply */
double* matrix_mult_vector(double* A, double* x_result, 
                    int row_count, int col_count) {
    int N = row_count; //assume equal to col_count as well
    double* C = (double*) malloc(sizeof(double) * N);
    double* rowA;
    double* rowB;
    double* rowC;
    double  sum = 0;
    int     i, j, k;
    for (i = 0; i < row_count; ++i) {
        rowA = &(A[i * N]);
        for (j = 0; j < 1; ++j) {
            sum = 0;
            rowB = &(x_result[j * N]);
            for (k = 0; k < col_count; ++k) {
                sum += rowA[k] * rowB[k]; 
            }
            C[i] = sum;
        }
    }
    return C;
}


void matrix_print(double* A, int row_size, int col_size) {
    printf("[ \n");
    for (int i = 0; i < col_size; ++i)
    {
        printf("    [ ");
        for (int j = 0; j <= row_size-2; ++j)
        {
            printf("%.2lf, ", A[ i * row_size + j ]);
        }
        printf("%.2lf ]", A[ i * row_size + row_size-1 ]);
        if (i != col_size-1) {
            printf(",\n");
        } else {
            printf("\n");
        }
    }
    printf("] \n\n");
}


void matrix_2d_print (double** B, int row_size, int col_size) {
    printf("[ \n");
    for (int i = 0; i < col_size; ++i)
    {
        printf("    [ ");
        for (int j = 0; j <= row_size-2; ++j)
        {
            printf("%.2lf, ", B[ i ][ j ]);
        }
        printf("%.2lf ]", B[ i ][ row_size-1 ]);
        if (i != col_size-1) {
            printf(",\n");
        } else {
            printf("\n");
        }
    }
    printf("] \n\n");
}





/* Some basic vector method definitons */

double* vector_new(int size) {
    if (size <= 0) {
        return NULL;
    }
    return (double*) calloc(size, sizeof(double));
}

void vector_init(double* vec, int vec_size, double value) {
    for (int i = 0; i < vec_size; ++i)
    {
        vec[i] = value;
    }
}

bool vector_delete(double* vec1) {
    if (vec1) {
        free(vec1);
    }
    return true;
}

double* vector_copy(double* vec1, int size) {
    double* vec2 = (double*) calloc(size, sizeof(double));
    for (int i = 0; i < size; ++i) {
         vec2[i] = vec1[i]; // [n] values set
    }
    return vec2;
}

/* vec1 - vec2.. 
    if retvec parameter was passed and it piitns to some allocated memory
        this this will set that memory instead.. 
    otherwise it creates a new vector on the heap..

    @return: returns a pointer to a new vector containing the results */
double* vector_subtract(double* vec1, double* vec2, 
                        int size,     double* retvec) {
    if (retvec == NULL) {
        retvec = (double*) calloc(size, sizeof(double));
    }
    for (int i = 0; i < size; ++i) {
        retvec[i] = vec1[i] - vec2[i];
    }
    return retvec;
}

double* vector_add( double* vec1, double* vec2, 
                    int size,     double* retvec) {
    if (retvec == NULL) {
        retvec = (double*) calloc(size, sizeof(double));
    }
    for (int i = 0; i < size; ++i) {
        retvec[i] = vec1[i] + vec2[i];
    }
    return retvec;
}

double* vector_scalar_mult( double* vec1, int size, 
                            double scalar_const, double* retvec) {
    if (retvec == NULL) {
        retvec = (double*) calloc(size, sizeof(double));
    }
    for (int i = 0; i < size; ++i) {
        retvec[i] = scalar_const * (vec1[i]);
    }
    return retvec;
}

double vector_dot(double* vec1, double* vec2, int size) {
    int product = 0;
    for (int i = 0; i < size; ++i) {
        product += vec1[i] * vec2[i];
    }
    return product;
}

/* given a function, draw onto the vector a set of values that correspond 
    to values in the range of function ffun on the domain i=0 to size
    ie. ffun should map element indecies into values.. */
bool vector_set_pattern( double* vec1, int size, 
                         double(*ffun)(double)) {
    if (!ffun) {
        return false;
    }
    for (int i = 0; i < size; ++i) {
        vec1[i] = ffun(i);
    }
    return true;
}

/* vec1 * M ..
    the matrix's column_size (number of rows) must be equal to vec_size
    (vector) times each column in the matrix */
double* vector_times_matrix (double* vec1,   int vec_size, 
                             double* matrix, int row_size) {
    double* mat = matrix_new(1, row_size);
    for (int i = 0; i < row_size; ++i) { // for all columns in rowsize number of columns
        mat[i] = vector_dot(vec1, matrix_col_get(matrix, i, row_size, vec_size), vec_size);
    }
    return mat;
}

void vector_print (double* vec1, int vec_size) {
    printf("[ ");
    for (int i = 0; i <= vec_size-2; ++i)
    {
        printf("%.2lf, ", vec1[i]);
    }
    printf("%.2lf ]\n", vec1[vec_size-1]);
}




void vector_randomize(double* b, int size){
    for (int i = 0; i < size; ++i) {
        b[i] = randomDouble();
    }
}