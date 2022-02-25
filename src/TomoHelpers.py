import numpy as np
from scipy.optimize import minimize 
import qutip
import pandas as pd



"""
get_state(fr)

Take a character ('H', 'V', 'D','A','R','L') and returns the vector corresponding to that basis state. E.g. input H outputs
np.matrix([[1], [0]])

Parameters:
----------------------
fr: string
    character corresponding to one of the qubit basis states. Available states are H, V, D, A, R, L.

    Note: Here, we use the convention where L represents the positive eigenstate of Pauli Y, and R represents the
    negative eigenstate.

"""

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

    else:
        raise ValueError("Not a valid basis state. Only valid inputs are 'H', 'V', 'D','A','R','L'.")
    
    return

"""
create_projections_str(n, curr_projs, base_strs)

Recursively produces len(curr_projs) ** 2 projections in string form. 

For example, if curr_projs =['H', 'V'], the function returns ['HH','HV','VH','VV']

Parameters:
----------------------
n: int
    number of qubits in the system.

curr_projs: list (str)
    the current list of projections (first input will always be an empty list)

base_strs: list (str)
    the list of strings from which you want to build the resulting list of projections (in the example, this would
    be ['H', 'V'])

"""
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

"""
str_to_state(str_projs, n)

Take a string composed of multiple basis states (e.g. "HH") as input, and returns the corresponding vector.
Parameters:
----------------------
str_projs: string
    string corresponding to one of the qubit basis states (any dimension). E.g. "HH" or "VR"

n: int
    number of qubits in the system.

"""

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


"""
create_t_matrix(n, params):

Consumes parammeters and returns the resulting T matrix for the system, upon which Rho is built. For more information on the 
structure of this matrix, check out J.B. Altepeter et al. (2019) 
(http://research.physics.illinois.edu/QI/Photonics/tomography-files/tomo_chapter_2004.pdf)


Parameters:
----------------------

n: int
    number of qubits in the system.

params: list (num)
    a list of length 4 ** n, representing the parameters that will put inserted into the matrix


"""

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

"""
maximum_liklihood_cost(params, n, counts, projections)

Computes and returns the current cost of our chosen cost function, which is defined in Altepeter et al. (2019). 


Parameters:
----------------------

params: list 
    paramaters used to compute our T matrix, and subsequently the rho matrix, which is used in the cost function.

n: int
    number of qubits in our system.

counts: list (num)
    List of counts for the tomography (list in the same order as the corresponding projections)

projections: list (np.matrix)
    List of the vector representing the states that were projected on for the tomography.

"""


def maximum_liklihood_cost(params, n, counts, projections):
    cost = 0 
    
    T = create_t_matrix(n,params)
    density_matrix = np.matmul(T.H, T)
    
    for i in range(len(projections)):
        
        predicted = np.real((((np.matmul(np.matmul(projections[i].H, density_matrix), projections[i]))).item(0)))
        cost += ((predicted - counts[i]) ** 2)/np.sqrt(predicted)
 
    return cost

"""
maximum_liklihood_estimation(n, counts, projections)

Performs the maximum likelihood technique on tomography data and returns a numpy array representing
the resulting density matrix.


Parameters:
----------------------

n: int
    number of qubits in our system.

counts: list (num)
    List of counts for the tomography (list in the same order as the corresponding projections)

projections: list (np.matrix)
    List of the vector representing the states that were projected on for the tomography.

"""

def maximum_liklihood_estimation(n, counts, projections):
    init_params = np.linspace(0.1, 1, 4 ** n)
    counts = counts/sum(counts)
    opt = {'disp':True,'maxiter':40000}

    soln_h = minimize(maximum_liklihood_cost, init_params, args=(n,counts,projections), method = 'SLSQP', options=opt)

    T = create_t_matrix(n,soln_h.x)
    density_matrix = 1/np.trace(np.matmul(T.H, T)) * np.matmul(T.H, T)

    return np.reshape(density_matrix, (2**n, 2**n))


"""
get_density_plots(rho,n)

Displays the 3D histogram of a density matrix

Parameters:
----------------------

rho:
    numpy array representing the density matrix of a quantum state

n:
    number of qubits in the system.

"""


def plot_density(n,rho):

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

"""
bloch_sphere(rho)

Displays the bloch sphere representation of a quantum state

Note: Only available for single qubit systems.

Parameters:
----------------------

rho:
    numpy array representing the density matrix of a quantum state

"""

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


"""
import_xl(n, file)

Reads in .xlsx file data and converts the information into counts and projections (to use this function, check out
example.xlsx to see how to format your data.)

Parameters:
----------------------

n: int
   number of qubits in systems

"""


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