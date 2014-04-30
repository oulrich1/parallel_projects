#include "approx_methods.h"


/* determines the minimum number of 
  trapezoids to use to approximate 
  the definite integral using approximate error.. */ 
approx_result_t
approx_min_interval(LD_T a, LD_T b, LD_T N, LD_T accepted_error){
    LD_T prev_approx  = integral_trap_method(&f1, a, b, N-1);
    LD_T approximate  = integral_trap_method(&f1, a, b, N);
    LD_T approx_error = abs_ld(approximate - prev_approx);
    while(approx_error > accepted_error){
        N++;
        prev_approx = approximate;
        approximate  = integral_trap_method(&f1, a, b, N);
        approx_error = abs_ld(approximate - prev_approx);
        debug("%.23Lf with N=%Lf, approx_error=%.14Lf\n",  approximate, N, approx_error);
    }
    debug("\n");
    approx_result_t p = {approximate, N, {approx_error, approx_error/approximate}};
    return p;
}

approx_result_t
approx_integral_find_min_traps_linearly(LD_T (*function_to_evaluate) (LD_T), 
                                    operation_constants_t   method_constants){
    LD_T a = method_constants.interval.a;
    LD_T b = method_constants.interval.b;
    LD_T accepted_error = method_constants.errors.relative;

    // mutable
    LD_T N = method_constants.interval.N;
    LD_T true_value  = method_constants.true_value;
    LD_T approximate = 0; 

    // the current relative error of the approx
    LD_T relative_error = 0; 

    approximate    = integral_trap_method(function_to_evaluate, a, b, N);
    relative_error = abs_ld(approximate - true_value)/abs_ld(true_value);   

    debug("--- Linear Search ------------------------------------------\n");
    debug("%.23Lf with N=%Lf, relative_error=%.20Lf\n",  approximate, N, relative_error);
    /* Just increment N by 1 until found.. */
    while(relative_error > accepted_error){
        N++;
        approximate  = integral_trap_method(function_to_evaluate, a, b, N);
        relative_error = abs_ld(approximate - true_value)/abs_ld(true_value);
        debug("%.23Lf with N=%Lf, relative_error=%.20Lf\n",  approximate, N, relative_error);
    }
    debug("\n");

    approx_result_t results = {approximate, N, {approximate - true_value, relative_error}};
    return results;
}

approx_result_t
approx_integral_find_min_traps_binary_method(LD_T (*function_to_evaluate) (LD_T), 
                                    operation_constants_t   method_constants){

    LD_T a = method_constants.interval.a;
    LD_T b = method_constants.interval.b;
    LD_T accepted_error = method_constants.errors.relative;

    // mutable
    LD_T N = method_constants.interval.N;
    LD_T true_value  = method_constants.true_value;
    LD_T approximate = 0; 

    // the current relative error of the approx
    LD_T relative_error = 0; 

    debug("--- Sub-interval Search ------------------------------------\n");

    // find min and max n to be used
    // as the interval for the exact find_n 
    LD_T left_n  = N-1;         // could just be 0.. but its the previous n by method of positive probing..
    LD_T right_n = N; // where the function is close enough.. but its an over estimation.. 
    LD_T factor = 2;  // speeds up linear probing.. by this factor

    approximate = 0;

    /* find the sub interval where there is probably a root */
    do{
        left_n = right_n;
        right_n *= factor;
        approximate = integral_trap_method(function_to_evaluate, a, b, right_n);
        debug("Left: %.15Lf.. Right: %.15Lf\n", left_n, right_n);
    } while(abs_ld(approximate - true_value)/true_value > accepted_error);
    
    N = left_n;

    debug("--- Binary Search ------------------------------------------\n");
        /* apply a binary search / midpoint method to narrow down the sub interval
        and find the root.. */
    operation_constants_t midpoint_constants = {{method_constants.interval.a, 
                                                 method_constants.interval.b}, 
                                                {0, method_constants.errors.relative}, 
                                                method_constants.true_value};

    N = approx_find_min_n_midpoint_rule(function_to_evaluate, 
                                        left_n, 
                                        right_n, 
                                        midpoint_constants);


    approximate = integral_trap_method(function_to_evaluate, a, b, N);
    relative_error = abs_ld(approximate - true_value)/abs_ld(true_value);
    
    approx_result_t results = {approximate, N, {approximate - true_value, relative_error}};
    return results;
}



/* Returns the N such that f(left)-f(right) is significantly small such
   that the next possible N in [left,right] is within the error bounds..
   Sandwich Theorem  */
LD_T
approx_find_min_n_midpoint_rule( LD_T (*func) (LD_T), LD_T left_n, LD_T right_n,
                                 operation_constants_t method_constants) {

    // integral of f over the interval [a,b] 
    LD_T a = method_constants.interval.a;    
    LD_T b = method_constants.interval.b;    
    LD_T N = 0; 

    // f over interval of N 
    LD_T left  = left_n;    
    LD_T right = right_n;   

    LD_T accepted_relative_error  = method_constants.errors.relative;
    LD_T true_value = method_constants.true_value;

    LD_T mid_n = 0;
    LD_T current_integral_f_mid_n = 0;
    LD_T current_error = 0;

    double prev_n_whole = 0;
    double n_whole      = 0;
    double n_residual   = 0;

    do{
        mid_n       = (right+left)/2; 
        current_integral_f_mid_n = integral_trap_method(func, a, b, mid_n);
        debug("left %.0Lf && right %.0Lf && int(f(x), mid_n) %.24Lf\n", 
                left, right, current_integral_f_mid_n);
        current_error = abs_ld(current_integral_f_mid_n - true_value)/true_value;
        if (current_error > accepted_relative_error) {
            left = mid_n;
            debug("Left  is now mid_n = %Lf\n", mid_n);
        } else if (current_error < accepted_relative_error) {
            right = mid_n;
            debug("Right  is now mid_n = %Lf\n", mid_n);
        } else{
            debug("Done...");
        }


        prev_n_whole = n_whole;
        n_residual = modf(mid_n, &n_whole); (void) n_residual;
        if (n_whole == prev_n_whole) {
            break; // DONE because the number is not changing significantly essntially.
        }
    } while (abs_ld( abs_ld(current_integral_f_mid_n- true_value)/true_value 
                    - accepted_relative_error) 
                    / accepted_relative_error
                > (accepted_relative_error)/10);

    N = mid_n;

    return N;
}



/* */
// same as above except the parameters are grouped together..
// also: this is based on the relative error instead..
approx_result_t
approx_integral_find_min_trapezoids(LD_T (*function_to_evaluate) (LD_T), 
                                    operation_constants_t   method_constants){

    approx_result_t results;
#ifdef APPROX_TRAP_COUNT_LINEAR
    results = approx_integral_find_min_traps_linearly(function_to_evaluate, method_constants);
#else 
    debug("--- Linear Search Disabled (re-enable in approx_methods.h) --- \n");
#endif



#ifdef APPROX_TRAP_COUNT_QUADRATIC_MIDPOINT
    results = approx_integral_find_min_traps_binary_method( function_to_evaluate, 
                                                            method_constants);

    method_constants.interval.N = results.N;

    operation_constants_t linear_constants = method_constants;

    results = approx_integral_find_min_traps_linearly(function_to_evaluate, linear_constants);
#else
    debug("--- Binary Search Disabled (re-enable in approx_methods.h) --- \n");
#endif

    return results;
}




/* Attempted to use bisection on a function defined
   to be the difference of relative errors.. :
   (abs(f(n) - true_value) / true_value) - accepted_relative_error = 0
   to find the root where this is true.. 
   however, the functions are not nicely mutable */
LD_T
approx_root_bisection(  LD_T (*function_to_evaluate) (LD_T), LD_T a, LD_T b, 
                        LD_T accepted_relative_error) {
    LD_T c = 0;
    LD_T f_c = 0;
    LD_T prev_f_c = 0;
    do
    {
        c = (b+a)/2;
        LD_T f_a = function_to_evaluate(a);
        LD_T f_b = function_to_evaluate(b);
        f_c = function_to_evaluate(c);
        if (f_a*f_c < 0) {  // one of the values is negative so there is a root
            b = c;
        } else if (f_c*f_b < 0) {  
            a = c;
        } else {   // neither are different signs.. 
            b = c; // THIS SHOULDNT HAPPEN unless we found a root..
        }
        if (abs_ld(prev_f_c - f_c) < accepted_relative_error 
            && abs_ld(f_c) > accepted_relative_error) {
            break;
        }
        prev_f_c = f_c;
    } while (abs_ld(f_c) > accepted_relative_error);
    if (abs_ld(f_c) > accepted_relative_error){
        debug("Bisection Warning: bad interval, f_c is not 0, bisection method broke..\n");
    }
    return c;
}



/* @params: report data type that holds all of the report data to write to a file.. */
void
approx_report_tuple_write_to_file(approx_report_tuple_t report_data, char* _filename){

    char  buffer[64];
    char  filename[64];  
    strcpy(filename, _filename);

    char  mode[2]   = "a+";
    FILE* ofp;
    ofp = file_create(filename, mode);

    // x1: procs
    //  x2: traps
    // x3: time
    sprintf(buffer, " (%Lf, %Lf), ", 
                report_data.x1, 
                //report_data.x2, 
                report_data.x3);

    file_write(ofp, buffer);

    // sprintf(buffer, "}");
    // file_write(ofp, buffer);

    file_kill(ofp);
}

