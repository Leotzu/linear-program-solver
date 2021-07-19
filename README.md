### CSC 445 - LP Project ###  
Author: Leo Mckee-Reid  

#### How to run the program: ####
* In terminal, write: python3 LP_solver.py < dir/file.txt
* Example: python3 LP_solver.py < test_LPs/input/vanderbei_exersise2.5.txt

#### Description of the program: #### 
* This is a linear program solver that uses the standard Simplex Method to find a solution.
* To choose a pivot, the program uses the largest coefficient method (choosing the largest coefficient out of all the objective function coefficients).
* If there is a negative value on the right hand side of one of the contraint equations, then the auxiliary method is used to find a feasible vertex before it begins making simplex pivots.
* If there is a chance of cycling (judged by a dictionary containing 2 or more zeros coefficients in the basis equations), then the lexicographic method is used to avoid cycling. This is done by assigning symbolic perturbations to each of the basis functions. These symbolic perturbations are used to choose the leaving variable in a strategic way, which ensures that a dictionary cannot cycle.

