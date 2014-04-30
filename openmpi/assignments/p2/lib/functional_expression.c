#include "functional_expression.h"

// hard coded functions
LD_T
f1_prime(LD_T x){
    return (-1/3)*sin(x/3) + (2/5)*sin(x/5) + (5/4)*cos(x/4);
}

LD_T
f1(LD_T x){
    return cos(x/3) - 2*cos(x/5) + 5*sin(x/4) + 8;
}

