import numpy as np
from scipy.optimize import minimize 
import qutip
import pandas as pd


def get_state(fr):
    if fr == "H" or fr == "h":
        return np.matrix([[1], [0]])
    if fr == "V" or fr == "v":
        return np.matrix([[0], [1]])
    if fr == "D" or fr == "d":
        return np.matrix([[1/np.sqrt(2)], [1/np.sqrt(2)]])
    if fr == "A" or fr == "a":
        return np.matrix([[1/np.sqrt(2)], [-1/np.sqrt(2)]])
    if fr == "R" or fr == "r":
        return np.matrix([[1/np.sqrt(2)], [-1j/np.sqrt(2)]])
    if fr == "L" or fr == "l":
        return np.matrix([[1/np.sqrt(2)], [1j/np.sqrt(2)]])
    
    return


def str_to_state(str_projs, n):
    
        projections = []
        for i in str_projs:

            if n == 1:
                state = get_state(i)
            for j in range(n - 1):

                curr_state = get_state(i[j])
                state = np.kron(curr_state, get_state(i[j+1]))

            projections.append(state)
        return projections


def create_t_matrix(n, params):
   # print(n)
    t_matrix = np.zeros((2 ** n, 2 ** n), dtype = "complex_")
    counter = 0
    idx = 0
    while(counter != 2 ** n):
        for i in range(2 ** n):
            for j in range(2 ** n):
                if i == j + counter:
                    if i == j:
                        t_matrix[i][j] = params[idx]
                        idx += 1
                    else:
                        
                        t_matrix[i][j] = (params[idx] + 1j * params[idx+1])
                        idx += 2
        counter += 1
    return np.asmatrix(t_matrix)

def maximum_liklihood_cost(params, n,counts, projections):
    cost = 0 
    
    T = create_t_matrix(n,params)
    density_matrix = np.matmul(T.H, T)
    
    for i in range(len(projections)):
        
        predicted = np.real((((np.matmul(np.matmul(projections[i].H, density_matrix), projections[i]))).item(0)))
        cost += ((predicted - counts[i]) ** 2)/np.sqrt(predicted)
 
    return cost

def maximum_liklihood_estimation(n, counts, projections):
    init_params = np.linspace(0.1, 1, 4 ** n)
    counts = counts/sum(counts)
    opt = {'disp':True,'maxiter':40000}

    soln_h = minimize(maximum_liklihood_cost, init_params, args=(n,counts,projections), method = 'SLSQP', options=opt)

    T = create_t_matrix(n,soln_h.x)
    density_matrix = 1/np.trace(np.matmul(T.H, T)) * np.matmul(T.H, T)

    return np.reshape(density_matrix, (2**n, 2**n))

def get_density_plots(rho,n):

    qutip_rho_real = qutip.Qobj(np.real(rho))
    qutip_rho_imag = qutip.Qobj(np.imag(rho))



    xlabels = create_projections_str(n, ['H', 'V'],['H', 'V'])
    ylabels = create_projections_str(n, ['H', 'V'], ['H', 'V'])

    fig, ax = qutip.matrix_histogram(qutip.Qobj(rho), xlabels, ylabels, limits=[-1,1],colorbar=True)

    ax.set_title("Real Rho", size=16)
    ax.view_init(azim=-55, elev=40)

    fig, ax = qutip.matrix_histogram(qutip_rho_imag, xlabels, ylabels, limits=[-1,1])
    ax.set_title("Imaginary Rho", size=16)

    ax.view_init(azim=-55, elev=40)
    return qutip.Qobj(rho)


def bloch_sphere(rho):
    b = qutip.Bloch()

    qutip_rho = qutip.Qobj(np.asarray(rho))

    b.xlabel[0] = "A"
    b.xlabel[1] = "D"

    b.zlabel[0] = "H"
    b.zlabel[1] = "V"

    b.ylabel[0] = "R"
    b.ylabel[1] = "L"

    b.add_states(qutip_rho)
    b.show()

def create_projections_str(n, curr_projs, base_strs):
    
    if n == 1:
        return curr_projs
    
    next_projections = []
    for i in curr_projs:
        for j in base_strs:
            curr_state = i + j
            next_projections.append(curr_state)
    n = n - 1;
    
    
    if n == 1:
        return next_projections
        
    else:
        return create_projections_str(n, next_projections, base_strs)


def plot_density(n, rho):
    qutip_rho_real = qutip.Qobj(np.real(rho))
    qutip_rho_imag = qutip.Qobj(np.imag(rho))

    xlabels = create_projections_str(n, ['H', 'V'],['H', 'V'])
    ylabels = create_projections_str(n, ['H', 'V'], ['H', 'V'])

    fig, ax = qutip.matrix_histogram(qutip.Qobj(rho), xlabels, ylabels, limits=[-1,1],colorbar=True)

    ax.set_title("Real Rho", size=16)
    ax.view_init(azim=-55, elev=40)

    fig, ax = qutip.matrix_histogram(qutip_rho_imag, xlabels, ylabels, limits=[-1,1])
    ax.set_title("Imaginary Rho", size=16)

    ax.view_init(azim=-55, elev=40)
    return qutip.Qobj(rho)

def import_xl(n, file):
    data = pd.read_excel(file, dtype = {"states": str, "counts":float})
    df = pd.DataFrame(data, columns= ['states', 'counts'])
    counts = np.zeros(len(df.values))
    projections_str = []

    for i in range(len(df.values)):
        projections_str.append(df.values[i][0])
        counts[i] = df.values[i][1]
    
    projections = str_to_state(projections_str, n)
    return projections, counts