import numpy as np
import time

# ========================= #
# Slack Inequal Constraints #
# ========================= #

def slack(A,b,c,hi):
    '''
    Input:
    A - (m*n ndarray) Assume nonsingular
    b - (1*m ndarray)
    c - (1*n ndarray) Coefficient of the cost function
    hi - (1*n ndarray) Upper Bounds
    '''
    (m,n) = A.shape
    for i in range(m):
        if b[i]<0:
            b[i] = -b[i]
            A[i,:] = -A[i,:]
    
    As = np.hstack((A,np.zeros((m,n))))
    Aslow = np.identity(n)
    Aslow = np.hstack((Aslow,np.identity(n)))
    As = np.vstack((As,Aslow))

    bb = b.reshape(m,1)
    hihi = hi.reshape(n,1)
    bs = np.vstack((bb,hihi))
    bs = list(bs.T.tolist()[0])

    cs = np.hstack((c,np.zeros(n)))
    return As,bs,cs,n


def revised_simplex(A,c,x,B,t):
    '''
    Input:
    A - (m*n ndarray) Assume nonsingular
    c - (1*n ndarray) coefficient of the cost function
    x - (1*n ndarray) a bfs
    B - (m list) the indice of the basis
    '''

    (m,n) = A.shape
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

        if enter == -1:             # Optimal
            optval = np.dot(c,x)
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
            status = -1             # Unbounded Below
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
        print("Iterations: ", t)

def PhaseI(A,b):
    '''
    Input:
    A - (m*n ndarray) Coefficient Matrix
    b - (m list) Constraints
    '''
    
    # Construct the Auxiliary Problem
    (m,n) = A.shape
    AP = A.copy()
    B = list(np.zeros(m)-1)
    for j in range(n):
        Ajnz = A[:,j].nonzero()[0]
        if len(Ajnz) == 1:
            Bi = Ajnz[0]
            B[Bi] = j
            if A[Bi,j] != 1:
                AP[Bi,:] = (1/A[Bi,j]) * A[Bi,:]
    
    basis = np.identity(m)
    for i in range(m):
        if B[i] == -1:
            AP = np.hstack((AP,basis[:,i].reshape(m,1)))
            B[i] = AP.shape[1]-1
        
    cp = np.hstack((np.zeros(n),np.ones(AP.shape[1]-n)))
    u = np.zeros(AP.shape[1])
    u[B] = b

    # First Simplex
    (status,iters,opv,xp,B,B_inv) = revised_simplex(AP,cp,u,B,0)

    # Check Feasibility
    if opv != 0:
        status = 0
        return status,iters,A,xp[:n],B
    
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
        print("Solution: x=",sol)
        print()
        return 1,opv,iters,rtime

    if status == -1:
        print("State: ", "Unbounded")
        print("Optimal Value =",opv)
        print("Total iterations =",iters)
        print("Running Time =",format(rtime, '.4f'),"s")
        print()
        return 0,opv,iters,rtime









