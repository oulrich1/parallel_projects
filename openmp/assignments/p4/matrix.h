#ifndef MATRIX_H
#define MATRIX_H 


#include <time.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
// #include <unistd.h> /* machine info */

#include <string.h>
#include <math.h>

#include <omp.h>


static const double MAX_VECTOR_VALUE = 1.0e10;

/*inline*/ double drand();
/*inline*/ double randomDouble();


/* Matrix Operations */

double* matrix_new      (int m, int n);

void    matrix_randomize (double* A, int row_size, int col_size);

bool    matrix_delete   (double* matrix);

bool    matrix_2d_delete(double** matrix2d);

double* matrix_copy     (double* A, int size);

/* copy 1d to a 2d.. */
double** matrix_copy_1d2d(double* A, int row_size, int col_size);

void    matrix_copy_columns (double* A, int row_size,   int col_size, 
                             double* B, int b_row_size );

bool    matrix_row_set  (double* A, int row_id, double* row_vector, int row_size);

bool    matrix_col_set  (double* A, int col_id,  int row_size, double* col_vector, int col_size);

double* matrix_row_get  (double* A, int row_id, int row_size);

double* matrix_col_get  (double* A, int col_id, int row_size, int col_size);

double* matrix2d_col_get(double** A, int col_id, int row_size, int col_size);

void    matrix_print    (double* A, int m, int n);

void matrix_2d_print (double** B, int m, int n);

double* matrix_augmented_new (double* A, double* b);


double* matrix_mult_vector(double* A, double* x_result, 
                            int row_count, int col_count);



/* Vector Operations */

double* vector_new          (int size);

void    vector_init         (double* vec, int vec_size, double value);

bool    vector_delete       (double* vec1);

double* vector_copy         (double* vec1, int size);

double* vector_subtract     (double* vec1, double* vec2, int size, double* retvec);

double* vector_add          (double* vec1, double* vec2, int size, double* retvec);

double* vector_scalar_mult  (double* vec1, int size, double scalar_const, double* retvec);

double  vector_dot          (double* vec1, double* vec2, int size);

bool    vector_set_pattern  (double* vec1, int size, double(*ffun)(double) );

double* vector_times_matrix (double* vec1,   int vec_size, double* matrix, int row_size);

void    vector_print        (double* vec1, int vec_size);

void    vector_randomize (double* b, int size);


#endif