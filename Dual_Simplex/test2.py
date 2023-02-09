# =================================================================== #
#                                                                     #
#                    Dual Simplex Method (Revised)                    #
#                                                                     #
# =================================================================== #

# [Small Examples] Comparison between Primal and Dual

'''
from Dual_Simplex_01 import *

A = np.array([[1,1,1,0],[1,0,0,1]])
b = np.array([[2],[1]])
c = np.array([-2,-1,0,0])
B = [1,4]

Dual_simplex(A,b,c,B)
'''

'''
from TwoPhase_Simplex_01 import *

A = np.array([[1,1,1,0],[1,0,0,1]])
b = np.array([2,1])
c = np.array([-2,-1,0,0])

TwoPhase_Simplex(A,b,c)
'''

'''
from Dual_Simplex_01 import *

A = np.array([[1,3,0,4,1],[1,2,0,-1,1],[-1,-4,3,0,0]])
b = np.array([[2],[2],[1]])
c = np.array([2,3,3,1,-2])
B = [3,4,5]
Dual_simplex(A,b,c,B)
'''

'''
from Dual_Simplex_01 import *

A = np.array([[1,2,3,0,1,0,0,0],[-1,2,6,0,0,1,0,0],[0,4,9,0,0,0,1,0],[0,0,3,1,0,0,0,1]])
b = np.array([[3],[2],[5],[1]])
c = np.array([0,0,0,0,1,1,1,1])
B = [5,6,7,8]
Dual_simplex(A,b,c,B)
'''

''' 
from Dual_Simplex_02 import *

A = np.array([[3,-1,-1,0],[1,0,0,-4],[-3,2,1,2]])
b = np.array([[-3],[-2],[6]])
c = np.array([-1,-2,-1,0])
B = [1,2,3]
Dual_simplex(A,b,c,B)
'''

''' 
from Dual_Simplex_02 import *

A = np.array([[3,-1,-1,0,1,0,0],[1,0,0,-4,0,1,0],[-3,2,1,2,0,0,1]])
b = np.array([[-3],[-2],[6]])
c = np.array([1,2,1,0,0,0,0])
B = [5,6,7]
Dual_simplex(A,b,c,B)
'''

# [Larger Scale] Comparison Between Primal and Dual
'''
from Dual_Simplex_01 import *

solvable = [15, 19, 26, 41, 98]
itersdual = []
for i in solvable:
    np.random.seed(i)
    m,n = 100,200
    A = np.random.randn(m,n)*10
    b = np.random.randn(m,1)*10
    c = np.random.randn(n)*10
    B = list(range(1,m+1))
    iter = Dual_simplex(A,b,c,B)
    itersdual.append(iter)

from TwoPhase_Simplex_01 import *

itersprimal = []
for i in solvable:
    np.random.seed(i)
    m,n = 100,200
    A = np.random.randn(m,n)*10
    b = np.random.randn(m)*10
    c = np.random.randn(n)*10
    state,opv,iters,rtime = TwoPhase_Simplex(A,b,c)
    itersprimal.append(iters)

print('---------------Comparison---------------')
print("Dual Simplex Method")
print("Average Iterations: ",np.mean(itersdual))
print("Primal Simplex Method")
print("Average Iterations: ",np.mean(itersprimal))
'''

# Comparison with Gurobi
'''
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp


solvable = [15, 19, 26, 41, 98]

for i in solvable:
    
    np.random.seed(i)
    m,n = 100,200

    A = np.random.randn(m,n)*10
    b = np.random.randn(m)*10
    c = np.random.randn(n)*10

    try:

        model = gp.Model("model")

        x = model.addMVar(shape=n, name="x")
        model.setObjective(c @ x, GRB.MINIMIZE)
        model.addConstr(A @ x == b, name="e1")

        loA = np.identity(n)
        lo = np.zeros(n)
        model.addConstr(loA @ x >= lo, name="ine1")

        # Optimize model
        model.optimize()

        print('Obj: %g' % model.ObjVal)

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))

    except AttributeError:
        print('Encountered an attribute error')
'''
