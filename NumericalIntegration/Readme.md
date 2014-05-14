Go into writeup.md to view the project design, description, my thoughts,
and implementation details. It also docuemnts the current TODOs and shortcomings..

The p2_report.pdf contains the required sumamry and results of the analysis of the
data that had been collected with this program.
--------------------------------------------------------------------------------

The integral methods and the functions themselves are coded in the lib/ 
The methods of approximation are also in the lib/ in "approx_methods.c"

The main program that implements mpi is within the test/ under p2/

See test/ for the main controllers..

--------------------------------------------------------------------------------


just run :

    "./run 10"

or run:

    "./run_test #procs [interval_a | interval_b | N_estimation]"

to run the parallel_integration on 10 processors, for example
If this doesnt work, continue reading..


--------------------------------------------------------------------------------


To build, run make:

    "make clean"
    "make"

This should build both programs.. but they should already be compiled

Then use 

    "mpirun -n $processor_count $target [...]"

To run the target program, either "find_min_trapezoids" or "parallel_integration"

