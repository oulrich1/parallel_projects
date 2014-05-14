#ifndef INTEGRAL_EXPRESSION_H
#define INTEGRAL_EXPRESSION_H

#include "util.h"

typedef struct 
{
    LD_T a;
    LD_T b;
    LD_T N;
    LD_T interval_width; 
} interval_t;

/* hardcoded functions */
LD_T f1_prime(LD_T x);
LD_T f1(LD_T x);

struct functional_t;

#endif