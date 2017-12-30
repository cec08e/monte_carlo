from __future__ import print_function
import numpy as np

# Ising grid nxm size
# enumerate states
#calculate partition
# mean magnetization per spin
#specific heat per son

class Ising2D(object):

    # Boltzmann constant in J/K
    b_const = 1.38e-23

    def __init__(self, n, m, ext_B = 0):
        ''' :param n: Size of first dimension of lattice.
            :type n: int.
            :param m: Size of second dimension of lattice.
            :type m: int.
            :param ext_B: Magnitude of external magnetic field.
            :type ext_B: float.
        '''
        self.n = n
        self.m = m
        self.lattice = [[-1 for i in range(self.m)] for j in range(self.n)]
        self.ext_B = ext_B

    def calc_Z(self, T):
        ''' Calculate partition function for lattice.
            :param T: Temperature.
            :type T: float.
        '''
