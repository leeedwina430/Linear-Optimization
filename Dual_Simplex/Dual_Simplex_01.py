# =================================================================== #
#                                                                     #
#                    Dual Simplex Method (Revised)                    #
#                           Gongxikuiyou                              #
#                                                                     #
# =================================================================== #
"""
    Expect Standard Form LP:
        minimize      c'x
        subject to   Ax = b
                     x >= 0
    
    Example:
        A = np.array([[1,2,3,0],[-1,2,6,0],[0,4,9,0],[0,0,3,1]])
        b = np.array([[3],[2],[5],[1]])
        c = np.array([1,1,1,0])
        Dual_simplex(A,b,c)
        
    Note:
        Dual_Simplex_01 : Revised Dual Simplex Method (Phase II)
        Dual_Simplex_02 : Version I with Devex Pricing. 
                          (Approximate Steepest Edge Rule)
"""

import numpy as np
import time


def Dual_simplex(A,b,c,B):
    '''
    Input:
    A - (m*n ndarray) Assume nonsingular
    b - (m*1 ndarray) RHS constraints
    c - (1*n ndarray) coefficient of the cost function
    B - (m list) the indice of the basis(not the real index)
    '''

    B = list(np.array(B)-1)
    (m,n) = A.shape
    B_inv = np.linalg.inv(A[:,B])
    N = list(set(np.arange(n))-set(B))
    w = np.matmul(c[B],B_inv)
    r = c[N] - np.matmul(w,A[:,N])

    t = 0
    while True:

        # Check Optimality
        beta = np.matmul(B_inv,b)
        pv = min(beta)
        if pv >= 0:
            x = np.zeros(n).reshape(n,1)
            x[B,:] = beta
            pval = np.dot(c,x)[0]
            print("Current x is Optimal. OPTVal = {:.2f}".format(pval))
            print("Total Iterations = ",t)
            return t
        p = np.where(beta==pv)[0][0]    # p Exits
        
        # Check Infeasibility
        u = B_inv[p,:]
        y = np.matmul(u,A[:,N])
        if min(y) >= 0 :
            print("Primal infeasible")
            return 0

        # Exit
        ny = np.where(y<0)[0]
        nroy = -r[ny]/y[ny]
        gamma = -min(nroy)
        q = np.where(nroy==-gamma)[0][0]     # q Enters
        q = ny[q]       #
        # Update Reduced Cost
        r = r - gamma*y
        r[q] = -gamma

        # Change the Basis
        d = np.matmul(B_inv,A[:,N[q]]).reshape(m,1)
        tempm = np.hstack((beta,B_inv,d))
        tempm[p,:] = tempm[p,:]*(1/d[p])
        for i in range(m):
            if i != p:
                tempm[i,:] = tempm[i,:] - tempm[p,:] * d[i]
        beta = tempm[:,0]
        B_inv = tempm[:,1:m+1]
        B[p],N[q] = N[q],B[p]

        t += 1
        print("Iterations: ", t)







