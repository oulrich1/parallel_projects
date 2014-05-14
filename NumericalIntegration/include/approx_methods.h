#ifndef APPROX_METHODS_H
#define APPROX_METHODS_H

// #include <gsl/gsl_sf_bessel.h>
// #include <gsl/gsl_math.h>

#include "integral_operations.h"
#include "file_io.h"

#ifndef APPROX_TRAP_COUNT_LINEAR
//#define APPROX_TRAP_COUNT_LINEAR 
#endif

#ifndef APPROX_TRAP_COUNT_QUADRATIC_MIDPOINT
#define APPROX_TRAP_COUNT_QUADRATIC_MIDPOINT 
#endif

void* (*approx_method)(LD_T);

// interval_t defined as a type within
// functional_expression.h

/* Errors */
typedef struct
{
    LD_T absolute;
    LD_T relative;
} errors_t;

// returned..
/* Container for the return values 
   of an approximation operation,
   returned when an approximation 
   occurs.. we usually want to know the result of the 
   approximation, and the number of iterations or even 
   guesses N.. within some error */
typedef struct
{
    LD_T      approximate;
    LD_T      N; // # of trapazoids, (aka: t)
    errors_t  error;
} approx_result_t;

// sent..
/* Container for the data required by 
   the method being utilized.. during 
   an algorithm such as root finding
   or trapazoidal integrals */
typedef struct
{
    interval_t interval;
    errors_t   errors;
    LD_T       true_value;
} operation_constants_t; 


typedef struct 
{
  LD_T x1;
  LD_T x2;
  LD_T x3;
} approx_report_tuple_t;

/* determines the minimum number of 
  trapezoids to use to approximate 
  the definite integral */ 
approx_result_t
approx_min_interval(LD_T a, LD_T b, LD_T N, LD_T accepted_error);

approx_result_t
approx_integral_find_min_trapezoids(LD_T (*function_to_evaluate) (LD_T), 
                                    operation_constants_t   method_constants);


LD_T
approx_find_min_n_midpoint_rule( LD_T (*func) (LD_T), LD_T left_n, LD_T right_n,
                                 operation_constants_t method_constants);

//   finds the root of f_prime within interval
LD_T 
approx_root_bisection(LD_T (*function_to_evaluate) (LD_T), 
                           LD_T a, LD_T b, LD_T accepted_error); 


// LD_T (*) (LD_T)
// approx_find_lagrange_poly( LD_T (*func) (LD_T), point_t* points);

void
approx_report_tuple_write_to_file(approx_report_tuple_t report_data, char* _filename);


#endif
