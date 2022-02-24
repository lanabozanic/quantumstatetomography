import numpy as np

def purity(rho):
    purity = np.trace(np.matmul(rho,rho))
    return np.round(np.real(purity), 3)

def s_param(n, rho):

    X = np.matrix([[0,1], [1,0]])
    Z = np.matrix([[1,0], [0,-1]])

    M1 = -1/np.sqrt(2)*(Z+X)
    M2 = -1/np.sqrt(2)*(Z-X)

    s = 0

    if n == 2:
        ops = [np.kron(X, M1),np.kron(X, M2),np.kron(Z, M1),np.kron(Z, M2)]
        for i in ops:
            expval = np.trace(np.matmul(rho,i))
            s += abs(expval)

        return np.round(s, 3)

    else:
        print("S-parameter: This is only available for two-qubit systems. S-param was not calculated")

def concurrence(rho):
    eigvals, eigvecs = np.linalg.eig(rho)
    eigvals = np.sort(eigvals)[::-1]
    con = eigvals[0]
    for i in eigvals:
        if i == eigvals[0]:
            continue
        else:
            con = con - i
    
    return max(0, np.round(np.real(con), 3))

def tangle(concurrence):
    return concurrence ** 2