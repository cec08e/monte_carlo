from __future__ import print_function
from random import choice, randint, random
from numpy import exp, power
import logging
from matplotlib import pyplot, colors
from functools import reduce
from scipy.integrate import simps
from scipy.optimize import curve_fit

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
        if T == 0:
            T = .00001
        self.accept_ratios = [exp(-(1.0/T)*4*self.J), exp(-(1.0/T)*8*self.J)]
        # Calculate initial magnetization and energy values
        self.mag = self.calc_magnetization()    # initial magnetization
        self.energy =  self.calc_energy() # initial energy
        self.mag_vals.append(self.mag/self.lat_size)
        self.energy_vals.append(self.energy/self.lat_size)
        logging.info('Starting simulation.')
        for i in range(max_steps):
            self.step()
            #if self.step_num in [0,100000,200000,400000,600000,1000000,2000000,4000000]:
            #    self.visualize_lattice()


    def calc_energy(self):
        # Calculate the instantaneous energy of the lattice
        # To sum over nearest neighbors and avoid double-counting,
        # we sum over all sites and only calc east and south neighbors.
        energy = 0
        for i in range(self.rows):
            for j in range(self.columns):
                energy += -self.J*self.lattice[i][j]*(self.lattice[i][(j+1)%self.columns] + self.lattice[(i+1)%self.rows][j])
        logging.info("Initial energy: " + str(energy))
        return energy




    def calc_magnetization(self):
        # Calculate the instantaneous magnetization of the lattice
        mag = 0
        for row in self.lattice:
            for spin in row:
                mag += spin
        logging.info("Initial magnetization: " + str(mag))
        return mag


    def step(self):
        #logging.info("Step " + str(self.step_num))
        # Select a new state by randomly choosing a spin to flip
        chosen_site_row = randint(0, self.rows - 1)
        chosen_site_col = randint(0, self.columns - 1)
        #logging.info("Site selected: " + str((chosen_site_row, chosen_site_col)))
        # Calculate the difference in energy between new and old state
        # Using the summation trick of Newman, Barkema (equation 3.10)
        delta_E = self.calc_delta_E(chosen_site_row, chosen_site_col)
        #logging.info("Delta E: " + str(delta_E))
        accept_flag = False
        if delta_E > 0:
            rand = random()
            if delta_E > self.min_positive:
                if rand > self.accept_ratios[1]:
                    # Accept the move with A = exp[-beta*delta_E]
                    #logging.info("Transition rejected.")
                    pass
                else:
                    #logging.info("Transition accepted.")
                    accept_flag = True
                    self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*-1
            elif rand > self.accept_ratios[0]:
                #logging.info("Transition rejected.")
                pass
            else:
                #logging.info("Transition accepted.")
                accept_flag = True
                self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*-1
        else:
            # Accept the move with A = 1
            #logging.info("Transition accepted.")
            accept_flag = True
            self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*-1
        # Update and record modified magnetization and energy values
        # Old values are stored in self.energy_vals[self.step_num-1] and self.mag_vals[self.step_num-1]
        if accept_flag:
            self.energy = self.energy + delta_E
            self.mag = self.mag + 2*self.lattice[chosen_site_row][chosen_site_col]


        if self.step_num%self.lat_size == 0:
            self.energy_vals.append(self.energy/self.lat_size) # record energy per site
            self.mag_vals.append(self.mag/self.lat_size)   # record magnetization per site

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

    def plot_mag_energy_per_site(self):
        #print("Energy vals: ", self.energy_vals)
        #print("Mag vals: ", self.mag_vals)

        pyplot.plot([i for i in range(len(self.energy_vals))],self.energy_vals, 'r') # plot energy per site
        pyplot.plot([i for i in range(len(self.mag_vals))],self.mag_vals, 'b') # plot magnetization per site
        pyplot.ylabel('Energy (red) and magnetization (blue) per site')
        pyplot.xlabel('Steps per site')
        pyplot.show()

    def autocorrelate(self, sweeps):
        # Calc and plot magnetization autocorrelation function as a function of sweeps
        # Average magnetization per site
        avg_mag = reduce((lambda x,y: x+y), self.mag_vals)/len(self.mag_vals)
        print("Average magnetization per site: ", avg_mag)
        avg_mag_sq = power(avg_mag, 2)

        chi_0 = self.autocorrelate_int(avg_mag_sq, 0) # Initial autocorrelation function chi(0)
        x = [i for i in range(sweeps)]
        y = [self.autocorrelate_int(avg_mag_sq, i, norm=chi_0) for i in range(sweeps)]
        pyplot.plot(x, y, 'b')
        pyplot.ylabel('Magnetization autocorrelation $\chi(t)$')
        pyplot.xlabel('Steps per site')

        # Add exponential fit
        #params, curve = curve_fit(lambda t,a,b: a*exp(b*t), x, y)
        #print("a: ", params[0], " and b: ", params[1])
        #pyplot.plot(x, [params[0]*exp(params[1]*i) for i in x], 'g')

        pyplot.show()

    def autocorrelate_int(self,avg_mag_sq, t, norm = 1):
        # Performs autocorrelation integral at sweep t
        # Can only integrate over t'=0 to t'=(self.step_num/self.lat_size)-t ?
        t_prime = [i for i in range((int(self.step_num/self.lat_size) - t))] # x samples
        y_samples = [((self.mag_vals[i]*self.mag_vals[i+t]) - avg_mag_sq)/norm for i in range((int(self.step_num/self.lat_size) - t))]
        return simps(y_samples, t_prime)



if __name__ == "__main__":
    # To do: add command line args
    logging.basicConfig(filename="ising2d.log", level=logging.INFO, format='%(message)s')
    lat = Ising2D(20, 20, init_T=0)
    # Exp 1
    #lat.print_lattice()
    #lat.visualize_lattice()
    lat.simulate(max_steps = 10000000, T=2.4)
    #lat.print_lattice()
    #lat.visualize_lattice()

    # Exp 2
    #lat = Ising2D(100, 100, init_T=100)
    #lat.simulate(max_steps=100000000, T=2.0)
    #lat.plot_mag_energy_per_site()
    lat.autocorrelate(1000)
