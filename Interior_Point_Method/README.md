# Interior Points Method

Implementation of Interior Points Method of Mehrotra's Predictor-Corrector Algorithm version in MATLAB. Suitable for Linear Problems of standard form:
```
min   c' * x
st.     A*x = b
         x >= 0
```

-----------------------------------------------------------------------------
# Quick Example

```Input
A = [1 2 2 1 0 0; 2 1 2 0 1 0; 2 2 1 0 0 1];
b = [20; 20; 20];
c = [-10 -12 -12 0 0 0];
IPM(A,b,c');
```

```Expected Output
iter  1: residual = 4.65e-01 obj = -1.531159e+02
iter  2: residual = 2.89e-02 obj = -1.344155e+02
iter  3: residual = 4.75e-05 obj = -1.359980e+02
iter  4: residual = 4.75e-08 obj = -1.360000e+02
----------------------------------------
[IPM] Status: Solved. Obj = -136.00
[m, n] = [3 6], Iterarion = 4, CPU = 0.046875
```

-----------------------------------------------------------------------------
# Files Describition

> test.m
> Some examples. To see the results, one can run each part in the file; to see different version of IPM, on can uncomment different script lines of "IPM()", "IPM_cho()", and "IPM_modcho()". 

> IPM.m
> Direct implementation of Predictor-Corrector Interior Point Method, reference from Numerical Optimization chapter 14, algorithm 14.3.

> IPM_cho.m
> Interior Point Method Using built-in function "cho()" for Standard Cholesky Decomposition in MATLAB.

> rsolvable.mat
>  An feasible and bounded example with 500*1000 dimensional equality constraint matrix A. Note that the first example show in "test.m" includes the data in this file.

> scrs8.mat, e226.mat, nug08.mat
>  These are given data from elearning. Note all of the selected problems' As have full row rank.

-----------------------------------------------------------------------------
# Function Parameters

```Inputs
%   A - (m*n) Equality Constraints' Coefficient (Should be Full Row Rank Matrix)
%   b - (m*1) vector; Equality Constraints' RHS
%   c - (n*1) vector; Coefficients of Object Funciton
```

```Returns
%   x - (n*1) vector; Primal Solution in Last Iteration
%   y - (m*1) vector; Dual Solution in Last Iteration
%   s - (n*1) vector; Dual Solution in Last Iteration
%   iter - scalar; Number of Total Iterations
%   rtime - scalar; Total Running Time
```
(Note: These functions will print the state while running.)



-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
# **References**

1. Video Resources:
   * [https://www.youtube.com/watch?v=rSN5UvP0weE] - CMU, Javier Pena, Lecture 16, Interior Point Methods(Part 1).
2. Text Resources:
   * Numerical Optimization, Chapter 14.