#ifndef INTEGRAL_OPERATIONS_H
#define INTEGRAL_OPERATIONS_H 

#include "functional_expression.h"


/* approx integral with trap method 
   over [a,b] in n sub-intervals */
LD_T
integral_trap_method(LD_T (*f)(LD_T), LD_T a, LD_T b, LD_T n);

#endif