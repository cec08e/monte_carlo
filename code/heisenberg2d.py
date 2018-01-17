from __future__ import print_function
from random import uniform
from numpy import power,sqrt

class Heisenberg2D(object):
    def __init__(self, rows, cols, init_T = 0, B = 0):
        self.initialize(rows, cols, init_T, B)

    def initialize(rows, cols, init_T, B):
        '''
        Initializes lattice with spin configuration
        determined by temperature.

        Does this make sense for Heisenberg? Since this is not a classical model...

        '''
        if init_T:
            # Generate random spin configuration
            self.lattice = [[ self.gen_random_spin() for j in range(cols)] for i in range(rows)]

        else:
            # Generate aligned spin configuration, based on B
            if B < 0:
                # Align spins along negative z-axis
                self.lattice = [[(0,0,-1) for j in range(cols)] for i in range(rows)]
            else:
                # Align spins along positive z-axis
                self.lattice = [[(0,0,1) for j in range(cols)] for i in range(rows)]


    def gen_random_spin(self):
        '''
        Generates a random spin orientation using Marsaglia's
        method. First, we uniformly generate two numbers
        x1 and x2 from [-1,1], throwing away the case where
        their collective magnitude is >= 1. From these,
        we generate the points x,y,z which are uniformly
        distributed on the sphere.
        '''
        x1 = uniform(-1,1)
        x2 = uniform(-1,1)
        mag_sq = (power(x1,2) + power(x2,2))
        while mag_sq >= 1:
            x1 = uniform(-1,1)
            x2 = uniform(-1,1)
            mag_sq = (power(x1,2) + power(x2,2))

        return (2*x1*sqrt(1-mag_sq), 2*x2*sqrt(1-mag_sq), 1-2*mag_sq)

    def perturb(self, spin):
        '''
        Perturb a spin provided as a tuple (x,y,z).
        - Generate random spin vector.
        - 
        '''
