#!/bin/bash

# this run_Test/sh script simply stores the test input
# and executes the program using the test input..
# the program is run by calling the "run" script which 
# runs the program..

 processor_count_list=(3 4 8 14 20 42 70 120 147 )
#processor_count_list=(2 3 9 15 21)
trapezoid_count_list=(1000000 3000000 4000000 4200000 4709890 4900000 5000000 9000000)

interval_a=100
interval_b=600
true_value="4003.7209001513268265929"
trapezoid_count_min="4709890"
accepted_accuracy="2.4976e-14"

# for proc in ${processor_count_list[*]}; do
#     for trap in ${trapezoid_count_list[*]}; do
#         num=`expr $proc` # - 1`
#         eval "./run $num $interval_a $interval_b $trap"
#     done
# done    

for proc in ${processor_count_list[*]}; do
    num=`expr $proc` # - 1`
    eval "./run $num $interval_a $interval_b $trapezoid_count_min"
done    
