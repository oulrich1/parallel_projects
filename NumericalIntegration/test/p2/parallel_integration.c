#include "approx_methods.h" /* Includes the util.h */
#include "file_io.h"

#include <mpi.h>            /* For MPI functions, etc */ 


// TODO: use these to express 
// math functions and operations
// #include <gsl/gsl_sf_bessel.h>
// #include <gsl/gsl_math.h>

/* just a bundle for mpi to pass around 
   so that the final processor knows about 
   the different data */
typedef struct
{
    LD_T value;
    LD_T dt;
} results_t;

/* 
 * Compile:    
 *    mpicc -g -Wall -std=c99 -o mpi_hello main.c
 * Usage:        
 *    mpirun  -n<number of processes> ./mpi_hello
 */

const int  MAX_STRING   = 64; // 3145728 9437184 4709890
const LD_T N_CONST      = 4709890.00; // Trapazoid count approximation 6.2e-15
const LD_T TRUE_VALUE   = 4003.7209001513268265929;

const LD_T accepted_relative_error = 2.4976e-14;

int main(int argc, char const *argv[])
{
  //char  greeting[MAX_STRING];  /* String storing message */
    int   comm_sz;               /* Number of processes    */
    int   my_rank;               /* My process rank        */

    LD_T N = N_CONST;

    LD_T x_0 = 100;
    LD_T x_n = 600;


    // LD_T root1 = approx_root_bisection(-100, 600, 0.001);
    // printf("x=%.23Lf, y=%.23Lf\n", root1, f(root1));

    /* Start up MPI */
    MPI_Init(NULL, NULL); 

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 

    /* Get my rank among all the processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

    if (comm_sz == 1) {
        exit(0);
    }

    if (my_rank == 0) { 

//#define INPUT_ASK_AB
#ifdef  INPUT_ASK_AB
        /* 0th process Master processor:
           perform user i/o */
        char c = 'Y';
        printf("You you like to enter [a,b] and N Trapezoids: [Y/n]\n");
        scanf("%s", &c);
        if (c != 'n') {
            printf("Please input [a,b]:\n");
            scanf("%Lf", &x_0);
            scanf("%Lf", &x_n);
            printf("Please input Trapezoid count (n):\n");
            scanf("%Lf", &N);
        }
#endif
    }

#ifndef  INPUT_ASK_AB
    if (argc >= 3) {
        x_0 = atoi(argv[2]);
    }
    if (argc >= 4) {
        x_n = atoi(argv[3]);
    }
    if (argc >= 5) {
        N   = atoi(argv[4]);
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);

    // N itself might not be seperatable into 
    // whole parts.. we need to find N whole parts
    //LD_T N_whole_traps = 0;

    LD_T dxPT = (x_n - x_0)/N;   // dx Per Trapazoid 

    // Traps per process = total traps / processes used to do work!
    LD_T T_per_process = N/(comm_sz-1);    
    

    LD_T dxPP  = T_per_process*dxPT;   // dx Per Processor


    if (my_rank != 0) { 
        debug("Process %d of %d, started!\n", my_rank, comm_sz);
        
        // process rank dependant
        // find the left bound of 
        // this processor's data
        LD_T a    = x_0 + (my_rank-1)*dxPP;         
        LD_T b    = a + dxPP;

        double time_start  = MPI_Wtime();
        LD_T   approximate = integral_trap_method(&f1, a, b, T_per_process);
        double time_end    = MPI_Wtime();

        const int SAMPLE_COUNT = 4;
        LD_T T1     = abs_ld(time_start - time_end);
        LD_T T_min  = T1;
        LD_T T_avg  = T1;
        for (int i = 0; i < SAMPLE_COUNT; ++i)
        {
            processor_results_t pr ={ my_rank, approximate, T_per_process, 0, 0, 
                                        abs_ld(time_end - time_start) };
            print_processor_results(pr);

            time_start = MPI_Wtime();
            LD_T approximate = integral_trap_method(&f1, a, b, T_per_process);
            approximate = approximate;

            time_end = MPI_Wtime();

            T1 = abs_ld(time_start - time_end);
            if(T1 < T_min){
                T_min = T1;
            }
            T_avg += T1;
        }
        T_avg = T_avg/SAMPLE_COUNT;

        printf("Rank#: %d : Finished (using t = %.0Lf"\
               "trapazoids): Minimum Time = %0.9Lf[sec]\n\n", 
                    my_rank, T_per_process, T_min);

        results_t results = {approximate, T_min};
        /* Send message to process 0 */ 
        MPI_Send(&results, sizeof(results_t), MPI_BYTE, 
                        0, 0, MPI_COMM_WORLD); 

        MPI_Barrier(MPI_COMM_WORLD);

    } else if (my_rank == 0) {  
        // TODO: parallelize this part as well...
        debug("Process %d of %d, Waiting for partial approximate solutions!\n", 
                my_rank, comm_sz);

        /* output the results to stdout and to file..  */
        MPI_Barrier(MPI_COMM_WORLD);

        results_t results = {0,0};
        LD_T total_approximate = 0;
        LD_T average_time = 0;
        for (int processor_id = 1; processor_id < comm_sz; processor_id++) {
            /* Receive message from process q */
            MPI_Recv(&results, sizeof(results_t), MPI_BYTE,  // MPI_BYTE
                    processor_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_approximate += results.value;
            average_time += results.dt; // temporarily store total dt
        } 
         
        // actually find average dt
        average_time = average_time / (comm_sz-1);

        LD_T relative_error = abs_ld(total_approximate-TRUE_VALUE)/TRUE_VALUE;

        printf("--- All Done -------------------------------------------\n");
        printf("Total Trapazoids:     %.10Lf\n", N);
        printf("Approx  Integral:     %.16Lf\n", total_approximate);
        printf("True    Integral:     %.16Lf\n", (LD_T)TRUE_VALUE);
        printf("Relative   Error:     %.16Lf\n", relative_error);

        printf("Average Elapsed Time: %.9Lf[sec]\n", average_time);
        printf("\n");

        if (relative_error <= accepted_relative_error) {
            printf( "Relative true error %0.16Lf is less than "\
                    "the acceptable error %0.16Lf\n", 
                    relative_error, accepted_relative_error);  
        } else {
            printf( "Relative true error %0.16Lf is NOT less than "\
                    "the acceptable error %0.16Lf\n", 
                    relative_error, accepted_relative_error);
        }

        char filename[64] = "points.out"; 
        approx_report_tuple_t report_tuple = {comm_sz-1, N, average_time};
        approx_report_tuple_write_to_file(report_tuple, filename);
    } 

    /* Shut down MPI */
    MPI_Finalize(); 

    return 0;
}

