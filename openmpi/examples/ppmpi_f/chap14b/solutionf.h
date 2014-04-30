C solutionf.h
C
C Declarations and definitions for manipulating the solutions
C
       integer SOLUTION_TAG  
       parameter (SOLUTION_TAG = 100)
       integer SolSizeMax
       parameter (SolSizeMax = 20) 
C
      integer  soln, rank
      parameter (soln = 0, rank = 1)
C
        integer best_solution(0:SolSizeMax)
       integer temp_solution(0:SolSizeMax)
C  first entry = cost, remaining entries, sibling   
C  ranks on path from root to solution node         
C  Thus solution_size = 1 + max_depth + 1           
       integer  solution_size
       common /Soln/ solution_size, 
     +      best_solution,temp_solution
C
C functions
       integer Best_solution_funct 
       integer Local_best_solution
       integer Initialize_soln
 
