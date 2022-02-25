import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from TomoHelpers import *
from RhoProperties import *


"""
class Rho()

a class defining the resulting density matrix of a tomography prodcedure, holding both the matrix and its key properites.

"""
class Rho():

    def __init__(self, n, rho):
        self.dims = [n,n]
        self.n = n
        self.matrix = rho
        self.purity = purity(self.matrix)
        if self.dims[0] == 2:
            self.s_param = s_param(self.n, self.matrix)
        self.concurrence = concurrence(self.matrix)
        self.tangle = tangle(self.concurrence)

    
"""
class Tomo()

object used to perform tomography

"""

class Tomo():

    def __init__(self, nqbits):
        self.n = nqbits


    """
    polarization_tomography_MLE()

    Performs the Maximum Likelihood Technique Algorithm to reconstruct the state. Returns a Rho() object.

    Parameters:
    ----------------------
    projections: list
        list of the states your tomography was projected on, in string form. e.g. ["H", "V", "D", "A", "L", "R"]
    
    counts: numpy array
        counts corresponding to the list of projections (e.g. if projections =[HH, HV...RR], your counts should be written in
        the same order)

    filename: string
        name of the the excel (.xlsx) file you wish to import your data from (to see how to format your data, check out
        example.xlsx)
    """

    def polarization_tomography_MLE(self, projections=[], counts=np.array([]), filename=None):

        if type(counts) != type(np.array([])):
            raise ValueError("Counts must be a numpy array")

        if filename != None:
            projections, counts = import_xl(self.n, filename)

        else:
            projections = str_to_state(projections, self.n)
        
        rho = maximum_liklihood_estimation(self.n, counts, projections)
        return Rho(self.n, rho)



    """
    bloch_visualization()

    Used to visualize a one-qubit state on the bloch sphere using QuTiP (only available for one-qubit tomography)

    Parameters:
    ----------------------

    rho: np.array
        matrix representing the density matrix of the state we wish to plot (must be 2x2)

    """

    def bloch_visualization(self, rho):
        if self.n != 1:
            raise ValueError("bloch_visualization is only available for 1 qubit systems")

        else:
            bloch_sphere(rho)


    """
    get_density_plot()

    Used to visualize any n qubit state on a histogram using QuTiP 

    Parameters:
    ----------------------

    rho: np.array
        matrix representing the density matrix of the state we wish to plot (must be 2**n x 2**n)

    """
    def density_plot(self, rho):
        if len(rho) != self.n * self.n:
            raise ValueError("rho must have length n * n. rho currently has length %d and n = %d" % (len(rho),self.n))
        plot_density(self.n, rho)

    """
    display_rho()

    Displays a formatted version of the rho matrix 

    Parameters:
    ----------------------

    rho: np.array
        matrix representing the density matrix of the state we wish to plot (must be 2**n x 2**n)

    """

    def display_rho(self,rho):
        return(qt.Qobj(rho))