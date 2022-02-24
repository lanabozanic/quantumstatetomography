import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from TomoHelpers import *
from RhoProperties import *

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

    

class Tomo():

    def __init__(self, nqbits):
        self.n = nqbits


    def polarization_tomography_MLE(self, projections=[], counts=[], filename=None):
        if filename != None:
            projections, counts = import_xl(self.n, filename)
        
        rho = maximum_liklihood_estimation(self.n, counts, projections)
        return Rho(self.n, rho)


    def bloch_visualization(self, rho):
        if self.n != 1:
            raise ValueError("bloch_visualization is only available for 1 qubit systems")

        else:
            bloch_sphere(rho)

    def get_density_plot(self, rho):
        if len(rho) != self.n * self.n:
            raise ValueError("rho must have length n * n. rho currently has length %d and n = %d" % (len(rho),self.n))
        plot_density(self.n, rho)

    def display_rho(self,rho):
        return(qt.Qobj(rho))