#include "util.h"

//##__VA_ARGS__
inline void debug_printf (const char* fmt, ...){
    if (!fmt){
        return;
    }

    va_list ap; /* arg pointer */
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
}

inline LD_T abs_ld(const LD_T val) { 
    return (val < 0) ? val * -1 : val; 
}

inline void
print_processor_results(processor_results_t processor_results) {
    #ifdef DEBUG_PROCESSOR_RESULTS
    printf("Rank#: %Lf:\n\
        approximate    =  %.16Lf;\n\
        t (# Traps)    =  %.2Lf;\n\
        true_error     =  %.20Lf;\n\
        relative_error =  %.20Lf;\n\
        run-time (sec) =  %.9Lf;\n\n", 
            processor_results.rank, 
            processor_results.approx, 
            processor_results.n, 
            processor_results.true_error, 
            processor_results.relative_error,
            processor_results.dt);
    #endif
}