Parallel/Numerical Writeup
--------------------------

Requirements:

    This assignment will do a couple of things in order to approximate the integral
    of a function. First it will find the number of trapezoids, N, required to 
    approximate the integral of the function within some error bounds. The size 
    of N will be directly proportional to the accuracy of the approximation, and 
    therefore inversely related to the magnitude of the error bounds. The next part
    of the "program" will use N to perform parallel approximations of the integral 
    and will also time these parallel computations and provide a report/summary of 
    computation.

Design:

    Use one program to calculate N in a serial program.
        The error bounds I looked at was the absolute value of the relative true 
        error. The goal was to find the N where the integral approximation was 
        equal to the true value for at least 14 significant digits.

        The true value of the integral will contain 20 significant digits.

        The N number of trapazoids will be found when the approximate integral is
        14 significant digits accurate to the true error. However, since I will be
        looking at the relative error of the approximation to the true value, I do
        not think this will guarentee accuracy within 14 sigificant digits since the
        base case for finding N will be when the relative error is less than the
        accepted relative error. 

        Code Example 1:
        accepted_error = 0.1e-14    // since we need 14 sig digs of the error to be 0
        until (relative_error <= accepted_error) {
            find N
        }

        Otherwise the base case would not depend on relative error and would simply 
        depend on the absolute error of the approx and true value. (absolute true error)

        Code Example 2:
        accepted_true_error = 0.1e(-14 + 4)    // 14 sig digits including 4 whole numbers
        until (abs(approx_value - true_value) < accepted_true_error) {
            find N
        }

        Howevever, this is not what I did. I just used the relative error as my 
        base case. I found the N where the relative error was less than or equal 
        to the accepted relative error.


    Now that I found my N number of trapazoids, i can find the integral in parallel
    using N as a constant.

    Using p = {1, 2, 4, 8, 10, 14, 100, 500, 1000} processes:
        Split up the intervalon which to calculate the integral into p sub intervals.
        Split N into p sub sets of n trapazoids. So now each sub interval has
            N/p number of trapazoids and an equal number of trapazoids for each subint.

        Time and Run the calculation in each process
        Repeat the calculattion at least 5 more times

        Display a report of the timings and approximations..
        The approximations should all be the same, although I think as more processors
        are used, there will be more possibility for round off and floating point 
        errors..

    After some sitting around, I realized it will be possible to quadratically
    probe the search space for an N. However this does not gaurentee a minimal 
    solution necessarily. So it would need to be a modified approach, where we 
    assume that under estimations and over estimations of N occur... However, 
    as soon as an over estimation occurs, our error will be accepted and therefore
    so will N. So, there are 2 clear ways to solve this issue:

        1.) Thinking about the curve fit of the set of points generated from the 
        integral estimation for each value of N.. as N -> inf .. the approximation
        gets asymtotically closer to the true value of the integral.
        What we want is N where the measured error is less than the accepted. 
        Error can only get better though.

        So, my idea is to curve fit the function of approximate values for some sample
        points of N.. then using that find the 0 of the function where 
        f(n) - error_bounds = 0... using a bisection method for example i could find
        the n = N where the approx is really close.. however i have to rethink this in terms
        of relative errors. So, the equation would be: 
            (abs(f(n) - true_value) / true_value) - accepted_relative_error = 0

        Curve Fitting ideas:
        
            could use some sort of curve fitting technique using genetic algorithm that 
            guesses polynomials over an interval.

            just use lagrange polynomials such that the error is.. some value of magnitude 
            proportional to the accepted_relative_error. P(n) + R(n) = f(n)

            R(n) is the error term that is necessary and it needs to be the absolute error..
            so R(n) = f(n) - true_value

        2.) Besides having to curve fit. I could just binary search until the sub 
        interval encapsulates a possible value for N where N is minimum and provides 
        enough accuracy. After finding the interval i could perform a binary search 
        in that interval to find the minimal value, or root, of the integral of the 
        function.. Since I am estimating with trapezoids this requires calculating 
        the integral for each value of N that I find. So runtime still has some overhead.

Program Implementaion:


    Integration of a function over an interval using N trapezoids.

    The program "find_min_trapezoids" finds the minimum number of trapezoids 
    necessary to approximate the integral within a determined relative error. 
    So that an approximation is accurate within 14 significant digits of the 
    actual true value, the relative error needs to be 1e-11 / true_value =~ 2.49e-15. 
    But i decided to use 5e-15 for the relative error instead. (I used 1e-11 as the 
        largest possible acceptable difference)

    After running find_min_trapezoids, i found N to be about 9437184.

    WHich is the value i hard coded into the next program "parallel_integration"
    The program "parallel_integration" integrations the function in parallel.


    DONE (TODO: (Issue: as n-processor count increases, the approximate is less accurate)
    Fix the error associated to splitting up the interval and trapezoid count
    into equal n-processor peices. There is error associated to chopping the 
    decimal point values in the trapaoid count (T_per_process) handed to each 
    processor. Process 0 needs to add the excess remainging whole number of 
    trapezoids and integrate over the last sub interval which contains the
    remaining number of trapezoids)
    
Current Issues:


    TODO: (Although the parallel work is super fast now.. 
        need to parallelize the summation..)
    As each processor complete's each share, need to sum the result in parallel. 
    Send the result to neighbor partner by rank who's value is power of two..
    from power = 1 to log n-processsors.
    

Todo:
    
    Graph the procedure of the "find_min_trapezoids" procesures.. 
        compare linear and quadratic times.
        compare the times for the approximation (should be linear)
    
    Graph the computation results of "parallel_integration" for each processor:
        such as time, approximation, N_p (# traps for that processor)

    Graph the final results of "parallel_integration" for a variable number of processors:
        for each number of processors in {1,2,4,6,8,10,50,80,120,200,250}
        display total averages:
            time, approximation, errors


