#ifndef UTIL_H
#define UTIL_H 

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* For strlen             */
#include <math.h>

#include <stdarg.h>

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

/* turtsmcgurts */

#define DEBUG_PROCESSOR_RESULTS
#define DEBUG_PRINTF

#ifdef DEBUG_PRINTF
    #define debug(fmt, ...) debug_printf(fmt, ##__VA_ARGS__)
#else
    #define debug(fmt, ...)
#endif


#ifndef uint
typedef unsigned int uint;
#endif

#ifndef LD_T
typedef long double LD_T;
#else
#warning "Typedef error:   cannot define 'long double' as 'LD_T'. LD_T already defined."
#warning "Typedef warning: defining LD_T as a double instead... could cause problems.."
typedef double LD_T
#endif

/* inline */ 
void debug_printf (const char* fmt, ...);

/* inline */ 
LD_T abs_ld(const LD_T val);

/* wrapper for the processor's results */
typedef struct 
{
    LD_T rank;
    LD_T approx;
    LD_T n;
    LD_T true_error;
    LD_T relative_error;
    LD_T dt;
} processor_results_t;

void
print_processor_results(processor_results_t processor_results);



#endif