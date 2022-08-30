import numpy as np
import time

def revised_simplex(A,c,x,B,t):
    '''
    Input:
    A - (m*n ndarray) Assume nonsingular
    c - (1*n ndarray) coefficient of the cost function
    x - (1*n ndarray) a bfs
    B - (m list) the indice of the basis
    '''
    #print("starting basis: B=",B)
    #print("starting bfs: x=",x)

    (m,n) = A.shape
    #print(A[:,B])
    B_inv = np.linalg.inv(A[:,B])
    
    while True:
        N = list(set(np.arange(n))-set(B))
        p = np.matmul(c[B], B_inv)
        c_ = c[N] - np.matmul(p, A[:,N])
        enter = -1
        mincj = np.infty

        # Find the enter index
        for j in range(len(N)):
            if c_[j] < 0:
                if c_[j] < mincj:
                    mincj = c_[j]
                    enter = N[j]

        if enter == -1:
            optval = np.dot(c,x)
            print("Current x is optimal.    x=",x)
            status = 1
            return status,t,optval,x,B,B_inv

        # Find the exit index
        u = np.matmul(B_inv, A[:,[enter]])
        l = -1
        theta = np.infty
        for i in range(m):
            if u[i] > 0:
                if (x[B[i]]/u[i]) < theta:
                    theta = x[B[i]]/u[i]
                    l = i
        
        if l == -1:
            #print("The problem is unbounded below.")
            status = -1
            return status,t,-np.inf,x,B,B_inv

        # Change the basis
        x[B[l]] = 0
        B[l] = enter
        B_inv[l,:] = B_inv[l,:]/u[l]
        x[B[l]] = theta
        for i in range(m):
            if i != l:
                B_inv[i,:] = B_inv[i,:] - u[i]*B_inv[l,:]
                x[B[i]] = x[B[i]] - theta*u[i]

        t += 1
        #print("Iterations: ", t, "    Current x=",x)
        print("Iterations: ", t)

def PhaseI(A,b):
    '''
    Input:
    A - (m*n ndarray) Coefficient Matrix
    b - (m list) Constraints
    '''

    # Construct the Auxiliary Problem
    (m,n) = A.shape
    for i in range(m):
        if b[i]<0:
            b[i] = -b[i]
            A[i,:] = -A[i,:]

    AP = np.hstack((A,np.identity(m)))
    u = np.hstack((np.zeros(n),b))
    B = list(range(n,n+m))
    cp = np.hstack((np.zeros(n),np.ones(m)))

    # First Simplex
    (status,iters,opv,xp,B,B_inv) = revised_simplex(AP,cp,u,B,0)

    # Check Feasibility
    if opv != 0:
        #print("The problem is infeasible.")
        status = 0
        return status,iters,A,xp,B
    
    # Drive Out the Artifitial Variables
    redundant = []
    for i in range(m):
        if B[i] >= n:
            for j in range(n):
                uu = np.matmul(B_inv,AP[:,j])
                if uu[i] != 0 and not (j in B):
                    B[i] = j
                    B_inv[i,:] = (1/uu[i])*B_inv[i,:]
                    B_invi = B_inv[i,:]
                    for k in range(m):
                        if k != i:
                            B_inv[k,:] = B_inv[k,:] - uu[k]*B_invi
                    break
            else:
                redundant.append(i)

    NA = np.zeros((1,n))
    NB = []
    for i in range(m):
        if i in redundant:
            continue
        else:
            NA = np.vstack((NA,A[i,:]))
            NB.append(B[i])

    x = xp[:n]
    status = 1
    return status,iters, NA[1:,:], x, NB

def TwoPhase_Simplex(A,b,c):
    start = time.perf_counter()

    print("Phase I")
    (status,t,NA,x,B) = PhaseI(A,b)
    
    rtime = time.perf_counter() - start
    if status == 0:
        print("-----------------------")
        print("State: ", "Infeasible")
        print("Optimal Value =",None)
        print("Total iterations =",t)
        print("Running Time =",format(rtime, '.6f'))
        #print("Solution: x=",x)
        print()
        return 0,None,t,rtime
    
    print("Phase II")
    (status,iters,opv,sol,Basis,B_inv) = revised_simplex(NA,c,x,B,t)

    rtime = time.perf_counter() - start
    print("-----------------------")
    if status == 1:
        print("State: ", "Solved")
        print("Optimal Value =", format(opv, '.6f'))
        print("Total iterations =",iters)
        print("Running Time =",format(rtime, '.6f'))
        #print("Solution: x=",sol)
        print()
        return 1,opv,iters,rtime

    if status == -1:
        print("State: ", "Unbounded")
        print("Optimal Value =",opv)
        print("Total iterations =",iters)
        print("Running Time =",format(rtime, '.6f'))
        #print("Solution: x=",sol)
        print()
        return 0,opv,iters,rtime









