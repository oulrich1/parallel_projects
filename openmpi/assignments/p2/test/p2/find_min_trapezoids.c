#include "approx_methods.h" /* Includes the util.h */

#include <mpi.h>            /* For MPI functions, etc */ 

// TODO: use these to express 
// math functions and operations
// #include <gsl/gsl_sf_bessel.h>
// #include <gsl/gsl_math.h>

/* 
 * Compile:    
 *    mpicc -g -Wall -std=c99 -o mpi_hello main.c
 * Usage:        
 *    mpirun  -n<number of processes> ./mpi_hello
 */

const LD_T T_approx     = 1; // Trapazoid count approximation
const LD_T TRUE_VALUE   = 4003.7209001513268265;

const LD_T accepted_relative_error = 2.4976e-14; 

int main(int argc, char const *argv[])
{
    LD_T x_0 = 100;
    LD_T x_n = 600;
    operation_constants_t operation_constants = {{x_0, x_n, T_approx}, {0, accepted_relative_error}, TRUE_VALUE};

    // playing with the bisection method on f1`
    // could remove this call, this isnt doing anything yet
    // printf("f`(%.16Lf) = 0\n\n", approx_root_bisection(&f1_prime, x_0, x_n, accepted_relative_error)); 

    double time_start = MPI_Wtime();
    // TODO: seperate the following functions 
    // into two functions: one that performs a linear searrch
    // and another that performs a quadratic + binary search
    approx_result_t solution = approx_integral_find_min_trapezoids(&f1, operation_constants);
    double time_end = MPI_Wtime();
        
    operation_constants.interval.N = solution.N;

    processor_results_t pr = {0, solution.approximate, 
                                solution.N,     
                                solution.error.absolute, 
                                solution.error.relative,
                                abs_ld(time_end - time_start) };

    print_processor_results( pr );

    return 0;
}


/* after running this program i get :  N = 9437184.00 */

// printf("Approximate = %0.16Lf.. Error = %0.16Lf.\n",  integral_trap_method(&f1, x_0, x_n, 1027008), 
//         (integral_trap_method(&f1, x_0, x_n, 1027008) - TRUE_VALUE) / TRUE_VALUE);
