from __future__ import print_function
from random import choice, randint
from numpy import exp


class Ising1D(object):

    def __init__(self, size, init_T = 0, B=0, J=1, bc='periodic'):
        self.lat_size = size
        self.bc = bc
        self.B = B
        self.J = J
        self.initialize(init_T)

        # Possible positive delta_E values are 4J only.
        # Therefore, our acceptance ratio for delta_E > 0
        # is always...
        self.accept_ratio = exp()

    def initialize(self, init_T = 0):
        # If T = 0, align all spins to ground state
        if init_T:
            # Assume T large enough to produce random spin state
            self.lattice = [choice([-1,1]) for i in range(self.lat_size)]
        else:
            if self.B < 0:
                # align all spins along -z
                self.lattice = [-1 for i in range(self.lat_size)]
            else:
                # align all spins along +z
                self.lattice = [1 for i in range(self.lat_size)]

    def simulate(self, T = 0):
        # Possible positive delta_E values are 4J only.
        # Therefore, our acceptance ratio for delta_E > 0
        # is always...
        self.accept_ratio = exp(-)

    def step(self):
        # Select a new state by randomly choosing a spin to flip
        chosen_site = randint(0, self.lat_size-1)
        # Calculate the difference in energy between new and old state
        # Using the summation trick of Newman, Barkema (equation 3.10)
        delta_E = 2*J*self.lattice[chosen_site]*(self.lattice[chosen_site+1] + self.lattice[chosen_site-1])
        if delta_E > 0:
            # Accept the move with A = exp[-beta*delta_E]
        else:
            # Accept the move with A = 1

    def print_lattice(self):
        for spin in self.lattice:
            if spin > 0:
                print('^', sep="", end="")
            else:
                print('v',sep="", end="")
        print("")

if __name__ == "__main__":
    lat = Ising1D(20, 9)
    lat.print_lattice()
