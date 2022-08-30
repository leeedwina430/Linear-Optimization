# Dual Simplex Method (Phase II)

## Interactive Files

> test.py
> It's a test file for the implementation. To see the results, one can delete ['''] at the head and bottom of each script. All the required libraries are well considered.

> TwoPhase_Simplex_0n.py
> These are two different versions of the Dual Simplex Method.
> Dual_Simplex_01 : Revised Dual Simplex Method (Phase II)
> Dual_Simplex_02 : Version I with Devex Pricing. 
>                      (Approximate Steepest Edge Rule)

## Users' Manual

### Expected problem

>  minimize  	     c'x
>
>  subject to 	    Ax = b
>
>  		    x >= 0

#### Reference

```python
from Dual_Simplex_01 import *
Dual_simplex(A,b,c,B)
# The answer will be printed 
```

#### Parameters:

    A - (m*n ndarray) Assume nonsingular
    b - (m*1 ndarray) RHS constraints
    c - (1*n ndarray) coefficient of the cost function
    B - (m list) the indice of the basis(not the real index)

#### Returns:

     iters - (int) The total number of iterations


#### Note:

This script will print the state while running.


#### Examples

```python
from Dual_Simplex_01 import *

b = np.array([[20],[20],[20]])
B = [3,4,5]
Dual_simplex(A,b,c,B)
```

```Output
Iterations:  1
Current x is Optimal. OPTVal = -3.00
Total Iterations =  1
```