from __future__ import print_function
import pickle
import shelve
from random import choice
from numpy import exp, power
from numpy.random import randint, rand, seed
import logging
from matplotlib import pyplot, colors
from functools import reduce
from scipy.integrate import simps
from scipy.optimize import curve_fit
import time

## TO DO: Separate simulation and analysis code.
## Add options for pickling data or using current run
## Store runs in increments of sweeps mapped to temperature
## E.g. {Temp1: [(mag1, e1),(mag2, e2),(mag3, e3)], Temp2: [(..),(..),(..)], ..}

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
        logging.info('Lattice size: ' + str(self.lat_size))
        logging.info('Initial T: ' + str(init_T))
        logging.info('B: ' + str(self.B))
        logging.info('J: ' + str(J))
        logging.info('Boundary conditions: ' + self.bc)
        self.initialize(init_T)
        #self.min_positive = 4*self.J


    def initialize(self, init_T = 0):
        self.sweep_num = 1 # initialize sweep count
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

    def simulate(self, num_sweeps, T, filename=None):
        self.mag_vals = []
        self.energy_vals = []
        #self.accept_ratios = [exp(-(1.0/T)*4*self.J), exp(-(1.0/T)*8*self.J)]

        if T == 0:
            T = .001
        # Calculate initial magnetization and energy values
        self.mag = self.calc_magnetization()    # initial magnetization
        self.energy =  self.calc_energy() # initial energy
        #logging.info("Initial magnetization: " + str(self.mag))
        #logging.info("Initial energy: " + str(self.energy))

        # mag_vals holds the magnetization per spin, calculated every sweep
        # energy_vals holds the total energy, calculated every sweep
        self.mag_vals.append(self.mag/self.lat_size)
        self.energy_vals.append(self.energy)
        #logging.info('Starting simulation.')
        for i in range(num_sweeps):
            self.sweep(T)
            if self.sweep_num%10000 == 0:
                print("Sweep # ", self.sweep_num)

        if filename:
            # Store the data associated with each sweep in a list of
            # tuples mapped to the temperature.
            data_file = shelve.open(filename)
            data_file[T] = list(zip(mag_vals, energy_vals))
            data_file.close()

    def calc_energy(self):
        # Calculate the instantaneous energy of the lattice
        # To sum over nearest neighbors and avoid double-counting,
        # we sum over all sites and only calc east and south neighbors.
        energy = 0
        for i in range(self.rows):
            for j in range(self.columns):
                energy += -self.J*self.lattice[i][j]*(self.lattice[i][(j+1)%self.columns] + self.lattice[(i+1)%self.rows][j])
        #logging.info("Initial energy: " + str(energy))
        return float(energy)

    def calc_magnetization(self):
        # Calculate the instantaneous magnetization of the lattice
        mag = 0
        for row in self.lattice:
            for spin in row:
                mag += spin
        #logging.info("Initial magnetization: " + str(mag))
        return float(mag)


    def sweep(self, T):
        # Perform as many steps as there are lattice sites
        for step in range(self.lat_size):
            # Select a new state by randomly choosing a spin to flip
            # Note: NumPy randint arguments are on a half-open interval
            chosen_site_row = randint(0, self.rows)
            chosen_site_col = randint(0, self.columns)

            # Calculate the difference in energy between new and old state
            # Using the summation trick of Newman, Barkema (equation 3.10)
            delta_E = self.calc_delta_E(chosen_site_row, chosen_site_col)

            #accept_flag = True
            # Double-check whether Python uses short circuit evaluation to
            # improve performance here
            if not ((delta_E > 0) and (rand() >= exp(-(1.0/T)*delta_E*self.J))):
                #accept_flag = False
                self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*(-1)
                self.energy = self.energy + delta_E
                self.mag = float(self.mag) + 2*float(self.lattice[chosen_site_row][chosen_site_col])
                #r = rand()
                #logging.info("random val is " + str(rand))
                #if delta_E > self.min_positive:
                #    if r < self.accept_ratios[1]:
                #        logging.info("Accepted. Rand is larger than larger exp.")
                #        accept_flag = True
                #        self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*(-1)
                #elif r < self.accept_ratios[0]:
                #    logging.info("Accepted. Rand is larger than smaller exp.")
                #    accept_flag = True
                #    self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*(-1)
            #else:
            #     Accept the move with A = 1
            #    logging.info("Accepted.")
            #    accept_flag = True
            # Update and record modified magnetization and energy values
            # Old values are stored in self.energy_vals[self.step_num-1] and self.mag_vals[self.step_num-1]
            #if accept_flag:
                #logging.info("Accept flag is true.")
                #self.lattice[chosen_site_row][chosen_site_col] = self.lattice[chosen_site_row][chosen_site_col]*(-1)
                #self.energy = self.energy + delta_E
                #self.mag = float(self.mag) + 2*float(self.lattice[chosen_site_row][chosen_site_col])
                #logging.info("New energy: " + str(self.energy))
                #print("new magnetization: ", self.mag, "             ", chosen_site_row, chosen_site_col)


            #if self.step_num%self.lat_size == 0:
            #    self.energy_vals.append(self.energy/self.lat_size) # record energy per site
            #    self.mag_vals.append(self.mag/self.lat_size)   # record magnetization per site
            #print(self.mag)

        # After sweep has completed, record energy and magnetization per spin
        self.energy_vals.append(self.energy)
        self.mag_vals.append(self.mag/self.lat_size)
        self.sweep_num += 1

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
        pyplot.imshow(self.lattice, cmap=cmap, norm=norm)
        pyplot.show()

    def plot_mag_energy_per_site(self):
        pyplot.plot([i for i in range(len(self.energy_vals))],[energy/self.lat_size for energy in self.energy_vals], 'r') # plot energy per site
        pyplot.plot([i for i in range(len(self.mag_vals))],self.mag_vals, 'b') # plot magnetization per site
        pyplot.ylabel('Energy (red) and magnetization (blue) per site')
        pyplot.xlabel('Steps per site (i.e. sweeps)')
        pyplot.show()

    def autocorrelate_mag(self, sweeps):
        # Calc and plot magnetization autocorrelation function as a function of sweeps
        # Average magnetization per site
        avg_mag = reduce((lambda x,y: x+y), self.mag_vals)/len(self.mag_vals)
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


    def autocorrelate_energy(self, sweeps):
        # Calc and plot energy autocorrelation function as a function of sweeps
        # Average energy per site
        avg_energy = reduce((lambda x,y: x+y), self.energy_vals)/len(self.energy_vals)
        avg_energy_sq = power(avg_energy, 2)

        chi_0 = self.autocorrelate_energy_int(avg_energy_sq, 0) # Initial autocorrelation function chi(0)
        x = [i for i in range(sweeps)]
        y = [self.autocorrelate_energy_int(avg_energy_sq, i, norm=chi_0) for i in range(sweeps)]
        pyplot.plot(x, y, 'b')
        pyplot.ylabel('Energy autocorrelation $\chi(t)$')
        pyplot.xlabel('Steps per site')

        # Add exponential fit
        #params, curve = curve_fit(lambda t,a,b: a*exp(b*t), x, y)
        #print("a: ", params[0], " and b: ", params[1])
        #pyplot.plot(x, [params[0]*exp(params[1]*i) for i in x], 'g')

        pyplot.show()

    def autocorrelate_energy_int(self, avg_energy_sq, t, norm = 1):
        # Performs autocorrelation integral at sweep t
        # Can only integrate over t'=0 to t'=(self.step_num/self.lat_size)-t ?
        #t_prime = [i for i in range((int(self.step_num/self.lat_size) - t))] # x samples
        t_prime = [i for i in range(self.sweep_num - t)] # x samples
        y_samples = [((self.energy_vals[i]*self.energy_vals[i+t]) - avg_energy_sq)/norm for i in t_prime]
        return simps(y_samples, t_prime)

    def autocorrelate_int(self,avg_mag_sq, t, norm = 1):
        # Performs autocorrelation integral at sweep t
        # Can only integrate over t'=0 to t'=(self.step_num/self.lat_size)-t ?
        #t_prime = [i for i in range((int(self.step_num/self.lat_size) - t))] # x samples
        t_prime = [i for i in range(self.sweep_num - t)] # x samples

        y_samples = [((self.mag_vals[i]*self.mag_vals[i+t]) - avg_mag_sq)/norm for i in t_prime]
        return simps(y_samples, t_prime)

    def suscept_v_temp(self, final_temp = 1, temp_step = .1, eq_time = None, cor_time = None):
        ''' Plots the susceptibility vs temperature for
            lattice. Plot begins with T=0 configuration and plots to final_temp
            in temp_step intervals.

            First, we calculate the susceptibility for the initial configuration.
            Then,
            - Simulate transition from T to T + temp_step (equilibration steps required).
            - Take measurements to find expectation value of magnetization per spin (correlation time required).
            - Calculate and plot the susceptibility.
            - Repeat steps above until we reach final_temp.

            For example, to calculate susceptibility of a 5x5 lattice
            from T=0 to T=5, we first create the lattice at temp 0.
            >>> lat = Ising2D(5, 5, init_T=0)
            Next, we call suscept_v_temp with the required parameters.
            >>> lat.suscept_v_temp(final_temp = 5, temp_step = .1, eq_time = 1250, cor_time = 1250)
            The suscept_v_temp function will handle all the required elibration and
            correlation simulation, but these step numbers must be supplied. For the
            5x5 example, a generous equilibration and correlation time is on the order
            of 50 sweeps (or 1250 steps). This can be determined qualitatively
            from the autocorrelation graph.

        '''

        # Begin by calculating and plotting the susceptibility for the initial
        # lattice configuration.
        # chi = (NT)(<m^2> - <m>^2)
        # Necessarily, our first plot point is 0
        chi_vals = [0]
        temp_vals = [0]
        curr_temp = 0

        while curr_temp < final_temp:
            curr_temp += temp_step
            self.simulate(eq_time, T = curr_temp)
            #mag_measurements = self.measure_mag_per_spin(cor_time, curr_temp)
            #print("mags: ", [item[0] for item in mag_measurements])
            self.simulate(cor_time, T = curr_temp)

            avg_mag_per_spin = reduce(lambda x,y: x+y, self.mag_vals)/len(self.mag_vals)
            #avg_mag_per_spin = reduce(lambda x,y: x+y, [item[0] for item in mag_measurements])/len(mag_measurements)
            #print("mags squared: ", [item[1] for item in mag_measurements])
            avg_mag_per_spin_sq = reduce(lambda x,y: x+y, [power(item,2) for item in self.mag_vals])/len(self.mag_vals)
            #print("<m>: ", avg_mag_per_spin)
            #print("<m^2>: ", avg_mag_per_spin_sq)
            #print("<m>^2: ", power(avg_mag_per_spin, 2))
            chi_vals.append((avg_mag_per_spin_sq - power(avg_mag_per_spin,2))*(self.lat_size*curr_temp))
            temp_vals.append(curr_temp)
            print("Current temp: ", curr_temp)

        pyplot.plot(temp_vals, chi_vals, 'bo')
        pyplot.ylabel('Susceptibility $\chi$')
        pyplot.xlabel('Temperature')
        pyplot.show()


    def spec_heat_v_temp(self, final_temp = 1, temp_step = .1, eq_time = None, cor_time = None):
        ''' Plots the specific heat per spin vs temperature for
            lattice. Plot begins with T=0 configuration and plots to final_temp
            in temp_step intervals.

            First, we calculate the specific heat for the initial configuration.
            Then,
            - Simulate transition from T to T + temp_step (equilibration steps required).
            - Take measurements to find expectation value of energy (correlation time required).
            - Calculate and plot the specific heat.
            - Repeat steps above until we reach final_temp.

            For example, to calculate specific heat per spin of a 5x5 lattice
            from T=0 to T=5, we first create the lattice at temp 0.
            >>> lat = Ising2D(5, 5, init_T=0)
            Next, we call spec_heat_v_temp with the required parameters.
            >>> lat.spec_heat_v_temp(final_temp = 5, temp_step = .1, eq_time = 1250, cor_time = 1250)
            The spec_heat_v_temp function will handle all the required elibration and
            correlation simulation, but these step numbers must be supplied. For the
            5x5 example, a generous equilibration and correlation time is on the order
            of 50 sweeps (or 1250 steps). This can be determined qualitatively
            from the autocorrelation graph. Specific heat will peak and change
            sharply because discrete energy shifts are large compared to overall energy?

        '''

        # Begin by calculating and plotting the specific heat for the initial
        # lattice configuration.
        # c = (1/NT^2)(<E^2> - <E>^2)
        # Necessarily, our first plot point is 0
        c_vals = [0]
        temp_vals = [0]
        curr_temp = 0

        while curr_temp < final_temp:
            curr_temp += temp_step
            self.simulate(eq_time, T = curr_temp)
            self.simulate(cor_time, T = curr_temp)

            #energy_measurements = self.measure_energy(cor_time, curr_temp)
            avg_energy = reduce(lambda x,y: x+y, self.energy_vals)/len(self.energy_vals)
            avg_energy_sq = reduce(lambda x,y: x+y, [power(item,2) for item in self.energy_vals])/len(self.energy_vals)
            c_vals.append((avg_energy_sq - power(avg_energy,2))/(self.lat_size*(power(curr_temp,2))))
            temp_vals.append(curr_temp)
            print("Current temp: ", curr_temp)

        pyplot.plot(temp_vals, c_vals, 'ro')
        pyplot.ylabel('Specific heat per spin $c$')
        pyplot.xlabel('Temperature')
        pyplot.show()


    def measure_energy(self, cor_time, curr_temp, num = 15):
        # Make energy measurements, num times
        energy_measurements = []

        for i in range(num):
            self.simulate(cor_time, T = curr_temp)
            energy = self.calc_energy()
            energy_measurements.append([energy, power(energy,2)])

        #print(energy_measurements)
        return energy_measurements

    def measure_mag_per_spin(self, cor_time, curr_temp, num = 30):
        # Make energy measurements, num times
        mag_measurements = []

        for i in range(num):
            self.simulate(cor_time, T = curr_temp)
            mag = self.calc_magnetization()  # when to divide by lat size?
            mag_measurements.append([mag/self.lat_size, power(mag,2)/self.lat_size])

        #print(energy_measurements)
        return mag_measurements

    # Analysis functions

    def mag_v_temp(self,  init_temp = 0, final_temp = 1, temp_step = .1, eq_time = None, cor_time = None, filename = None):
        ''' Plots the mean magnetization per spin vs temperature for
            lattice. Plot begins with T=0 configuration and plots to final_temp
            in temp_step intervals.
        '''
        #mag_vals = [1]
        mag_vals = []
        temp_vals = []
        curr_temp = init_temp

        # TO DO: Add support for pulling data from file.
        while curr_temp < final_temp:
            #self.initialize(init_T=100)

            #seed(int(time.time()))
            curr_temp += temp_step
            ##logging.info("*****************************")
            #logging.info("Now trying to come to equilibrium at T=" + str(curr_temp))
            #if curr_temp > 1 and curr_temp < 3:
            #    self.simulate(eq_time*10, T = curr_temp, calc_mag = False, calc_E = False)
            #else:
            self.simulate(eq_time, T = curr_temp)


            #mag_measurements = self.measure_mag_per_spin(cor_time, curr_temp)
            #avg_mag_per_spin = reduce(lambda x,y: x+y, [item[0] for item in mag_measurements])/len(mag_measurements)
            #if curr_temp > 1 and curr_temp < 3:
            #    if curr_temp > 1.5 and curr_temp < 2.5:
            #        self.simulate(cor_time*4, T = curr_temp, calc_mag = False, calc_E = False)
            #    else:
            #        self.simulate(cor_time*2, T = curr_temp, calc_mag = False, calc_E = False)

            #else:
            self.simulate(cor_time, T = curr_temp)

            avg_mag_per_spin = reduce(lambda x,y: x+y, self.mag_vals)/len(self.mag_vals)


            #print("Avg mag at temp ", curr_temp, ": ", avg_mag_per_spin)
            #print("mags squared: ", [item[1] for item in mag_measurements])
            #avg_mag_per_spin_sq = reduce(lambda x,y: x+y, [item[1] for item in mag_measurements])/len(mag_measurements)
            #print("<m>: ", avg_mag_per_spin)
            #print("<m^2>: ", avg_mag_per_spin_sq)
            #print("<m>^2: ", power(avg_mag_per_spin, 2))
            #chi_vals.append((avg_mag_per_spin_sq - power(avg_mag_per_spin,2))*(self.lat_size*curr_temp))
            mag_vals.append(avg_mag_per_spin)
            temp_vals.append(curr_temp)
            print("Current temp: ", curr_temp)

        pyplot.plot(temp_vals, mag_vals, 'go')
        pyplot.ylabel('Magnetization per spin $m$')
        pyplot.xlabel('Temperature')
        pyplot.show()





if __name__ == "__main__":
    # To do: add command line args
    logging.basicConfig(filename="ising2d.log", filemode='w',level=logging.WARNING, format='%(message)s')
    #lat = Ising2D(20, 20, init_T=0)
    # Exp 1
    #lat.print_lattice()
    #lat.visualize_lattice()
    #lat.simulate(max_steps = 10000000, T=2.4)
    #lat.print_lattice()
    #lat.visualize_lattice()

    # Exp 2
    #lat = Ising2D(100, 100, init_T=100)
    #lat.simulate(max_steps=100000000, T=2.0)
    #lat.plot_mag_energy_per_site()
    #lat.autocorrelate(1000)


    lat = Ising2D(50, 50, init_T=0)
    #lat.simulate(num_sweeps=10000, T=5)
    #lat.plot_mag_energy_per_site()
    #lat.autocorrelate_mag(800)
    #lat.autocorrelate_energy(500)

    #lat.spec_heat_v_temp(final_temp = 5, temp_step = .1, eq_time = 50000, cor_time = 600000)
    #lat.suscept_v_temp(final_temp = 10, temp_step = .1, eq_time = 50000, cor_time = 600000)
    ##lat.spec_heat_v_temp(final_temp = 5, temp_step = .1, eq_time = 300000, cor_time = 150000)

    ### 100x100 - eq_time = 1500000 cor_time = 150

    #lat.mag_v_temp(init_temp = 0, final_temp=5, temp_step=.1, eq_time=1000, cor_time=1000)
    lat.suscept_v_temp(final_temp=5, temp_step=.1, eq_time=1000, cor_time=1000)
