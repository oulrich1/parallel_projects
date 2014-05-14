#include "integral_operations.h"

/* approx integral with trap method 
   over [a,b] in n sub-intervals */
LD_T
integral_trap_method(LD_T (*f)(LD_T), LD_T a, LD_T b, LD_T n){

    // keep track of residual/fractional trapezoidal area..
    double n_exact = n;
    double n_whole = n;
    double n_residual = modf(n_exact, &n_whole); // of n_exact, the fractional bit is returned 
    (void)n_residual;

    LD_T dh = (b-a)/n;
    LD_T c  = a + n_whole*dh;

    LD_T f_a = f(a);
    LD_T f_c = f(c);
    LD_T f_b = f(b);

    LD_T x = a;
    
    // whole parts of the integral
    LD_T approx_sum_int = (f_a + f_c)/2;

    // perform the integral from 
    for (int i = 1; i <= n_whole-1; i++) {
        x = a + i*dh;
        approx_sum_int += f(x);
    }
    approx_sum_int *= dh;

    // now residual bit..
    x = c;
    approx_sum_int += ((f_c + f_b)/2)*n_residual*dh;

    return approx_sum_int;
}