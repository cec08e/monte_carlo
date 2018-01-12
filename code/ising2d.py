from __future__ import print_function
from random import choice, randint, random
from numpy import exp
import logging


class Ising2D(object):

    def __init__(self, rows, columns, init_T = 0, B=0, J=1, bc='periodic'):
        self.rows = rows
        self.columns = columns
        self.lat_size = rows*columns
        self.bc = bc
        self.B = B
        self.J = J
        self.step_num = 1
        logging.info('**************************')
        logging.info('Creating a 2D lattice with ' + str(self.lat_size) + ' sites.')
        logging.info('Rows: ' + str(self.rows))
        logging.info('Columns: ' + str(self.columns))
        logging.info('Initial T: ' + str(init_T))
        logging.info('B: ' + str(self.B))
        logging.info('J: ' + str(J))
        logging.info('Boundary conditions: ' + self.bc)
        self.initialize(init_T)


    def initialize(self, init_T = 0):
        # If T = 0, align all spins to ground state
        if init_T:
            # Assume T large enough to produce random spin state
            self.lattice = [[choice([-1,1]) for i in range(self.columns)] for j in range(self.rows)]
        else:
            if self.B < 0:
                # align all spins along -z
                self.lattice = [[-1 for i in range(self.columns)] for j in range(self.rows)]
            else:
                # align all spins along +z
                self.lattice = [[1 for i in range(self.columns)] for j in range(self.rows)]

    def simulate(self, T = 0):
        # Possible positive delta_E values are 2J and 4J only.
        # Therefore, our possible acceptance ratios for delta_E > 0 are
        self.min_positive = 2*J
        self.accept_ratios = [exp(-(1.0/T)*4*self.J), exp(-(1.0/T)*8*self.J)]
        logging.info('Starting simulation.')
        for i in range(1000):
            self.step()

    def step(self):
        logging.info("Step " + str(self.step_num))
        # Select a new state by randomly choosing a spin to flip
        chosen_site_row = randint(0, self.rows - 1)
        chosen_site_col = randint(0, self.colums - 1)
        logging.info("Site selected: " + str((chosen_site_row, chosen_site_col)))
        # Calculate the difference in energy between new and old state
        # Using the summation trick of Newman, Barkema (equation 3.10)
        delta_E = 2*self.J*self.lattice[chosen_site]*(self.lattice[(chosen_site+1)%self.lat_size] + self.lattice[(chosen_site-1)%self.lat_size])
        logging.info("Delta E: " + str(delta_E))
        if (delta_E > 0) and (random() > self.accept_ratio):
            # Accept the move with A = exp[-beta*delta_E]
            logging.info("Transition rejected.")
        else:
            # Accept the move with A = 1
            logging.info("Transition accepted.")
            self.lattice[chosen_site] = self.lattice[chosen_site]*-1
        self.step_num += 1


    def print_lattice(self):
        for spin in self.lattice:
            if spin > 0:
                print('^', sep="", end="")
            else:
                print('v',sep="", end="")
        print("")

if __name__ == "__main__":
    logging.basicConfig(filename="ising2d.log", level=logging.INFO, format='%(message)s')
    lat = Ising1D(20, init_T=0)
    lat.print_lattice()
    lat.simulate(T=2.4)
    lat.print_lattice()
