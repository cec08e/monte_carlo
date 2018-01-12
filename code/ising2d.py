from __future__ import print_function
from random import choice, randint, random
from numpy import exp
import logging
from matplotlib import pyplot, colors


class Ising2D(object):

    def __init__(self, rows, columns, init_T = 0, B=0, J=1, bc='periodic'):
        self.rows = rows
        self.columns = columns
        self.lat_size = rows*columns
        self.bc = bc
        self.B = B
        self.J = J
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
        self.step_num = 1 # initialize step count
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
        self.mag_vals = []
        self.energy_vals = []

    def simulate(self, max_steps, T = 0, calc_mag = True, calc_E = True):
        # If calc_mag or calc_E are true, magnetization and energy are calculated
        # and recorded over every sweep interval.
        # Possible positive delta_E values are 4J and 8J only.
        # Therefore, our possible acceptance ratios for delta_E > 0 are
        self.min_positive = 4*self.J
        self.accept_ratios = [exp(-(1.0/T)*4*self.J), exp(-(1.0/T)*8*self.J)]
        #mag =    # initial magnetization
        #energy =  # initial energy
        logging.info('Starting simulation.')
        for i in range(max_steps):
            self.step()
            #if self.step_num in [0,100000,200000,400000,600000,1000000,2000000,4000000]:
            #    self.visualize_lattice()
            #if self.num_steps % self.lat_size == 0:
                # We've completed another sweep
                # Very costly check - change condition
                #TO DO
                #pass


    def step(self):
        logging.info("Step " + str(self.step_num))
        # Select a new state by randomly choosing a spin to flip
        chosen_site_row = randint(0, self.rows - 1)
        chosen_site_col = randint(0, self.columns - 1)
        logging.info("Site selected: " + str((chosen_site_row, chosen_site_col)))
        # Calculate the difference in energy between new and old state
        # Using the summation trick of Newman, Barkema (equation 3.10)
        delta_E = self.calc_delta_E(chosen_site_row, chosen_site_col)
        logging.info("Delta E: " + str(delta_E))
        if delta_E > 0:
            rand = random()
            if delta_E > self.min_positive:
                if rand > self.accept_ratios[1]:
                    # Accept the move with A = exp[-beta*delta_E]
                    logging.info("Transition rejected.")
                else:
                    logging.info("Transition accepted.")
                    self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*-1
            elif rand > self.accept_ratios[0]:
                logging.info("Transition rejected.")
            else:
                logging.info("Transition accepted.")
                self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*-1
        else:
            # Accept the move with A = 1
            logging.info("Transition accepted.")
            self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*-1
        self.step_num += 1



    def calc_delta_E(self, row, col):
        # neighbor_sum is sum of nearest neighbor spin values
        neighbor_sum = self.lattice[(row-1)%self.rows][col] # north neighbor
        neighbor_sum += self.lattice[(row+1)%self.rows][col] # south neighbor
        neighbor_sum += self.lattice[row][(col-1)%self.columns] # west neighbor
        neighbor_sum += self.lattice[row][(col+1)%self.columns] # east neighbor
        return 2*self.J*self.lattice[row][col]*neighbor_sum


    def print_lattice(self):
        for row in self.lattice:
            for spin in row:
                if spin > 0:
                    print('^', sep="", end="")
                else:
                    print('v',sep="", end="")
            print("")

    def visualize_lattice(self):
        norm = colors.Normalize(vmin=-1, vmax=1)
        cmap = colors.ListedColormap(['white','black'])
        #bounds=[-2,0,.1,2]
        #norm = colors.BoundaryNorm(bounds, cmap.N)

        pyplot.imshow(self.lattice, cmap=cmap, norm=norm)
        pyplot.show()

if __name__ == "__main__":
    logging.basicConfig(filename="ising2d.log", level=logging.INFO, format='%(message)s')
    lat = Ising2D(100, 100, init_T=0)
    #lat.print_lattice()
    lat.visualize_lattice()
    lat.simulate(max_steps = 100000, T=2.4)
    #lat.print_lattice()
    lat.visualize_lattice()
