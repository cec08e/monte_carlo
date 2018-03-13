""" Cython Optimization """
# cython: profile=True

from __future__ import print_function
from numpy.random import uniform
from numpy.random import random, randint, rand
from numpy import power,sqrt,pi,sin,cos,arccos,arctan
from numpy import dot, cross, add, exp, linspace, array
from numpy.linalg import norm
from numpy import float64
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors
import logging
import json
import time
import pickle

# Change last dimension of lattice to 3 and all spin vectors to dimension 3?

cdef class Heisenberg3D(object):
    # Cython typing class attributes
    cdef int rows, cols, sweep_num
    cdef float J_intra, J_inter
    cdef float B, k1, k2, init_T
    cdef int lat_size
    cdef double [:,:,:,:] lattice

    def __init__(self, int rows, int cols, float init_T = 5, float B = 0, float J_intra = 1, float J_inter = .1, float k1 = 0, float k2 = 0):
        self.rows = rows
        self.cols = cols
        self.J_intra = J_intra     # Intra-layer interaction
        self.J_inter = J_inter     # Inter-layer interaction
        self.B = B
        self.k1 = k1
        self.k2 = k2
        self.lat_size = rows*cols
        self.init_T = init_T
        self.initialize(init_T, B)
        self.sweep_num=0

    def initialize(self, init_T, B):
        '''
        Initializes lattice with spin configuration
        determined by temperature.
        '''
        # Generate random spin configuration
        # Two layers
        self.lattice = array([[[ self.gen_random_spin() for j in range(self.cols)] for i in range(self.rows)] for g in range(2)])


    cpdef double[:] gen_random_spin(self):
        '''
        Generates a random spin orientation using Marsaglia's
        method. First, we uniformly generate two numbers
        x1 and x2 from [-1,1], throwing away the case where
        their collective magnitude is >= 1. From these,
        we generate the points x,y,z which are uniformly
        distributed on the sphere.
        '''
        cdef double x1
        cdef double x2
        cdef double mag_sq

        x1 = uniform(-1,1)
        x2 = uniform(-1,1)
        mag_sq = power(x1,2) + power(x2,2)
        while mag_sq >= 1:
            x1 = uniform(-1,1)
            x2 = uniform(-1,1)
            mag_sq = power(x1,2) + power(x2,2)

        return array([2*x1*sqrt(1-mag_sq), 2*x2*sqrt(1-mag_sq), 1-2*mag_sq], float)

    cdef double[:] perturb(self, double [:] spin):
        '''
        - Generate point x,y in circle with radius R
          see: http://xdpixel.com/random-points-in-a-circle/
        - Place point on tangent plane defined by current spin vector
            - Shift x,y origin by adding current spin vector.
            - Rotate plane by theta of current spin vector.

        '''
        cdef double R, r, arg, theta, phi
        cdef double x, y, z, u_x, u_y

        R = .6    # Optimal value should be about .4 according to Nehme et al.
        r = sqrt(random())*R
        arg = random()*2*pi   # can make 2pi a const

        # Express random point as cartesian coordinate in same reference frame as spin
        x = r*sin(arg)
        y = r*cos(arg)
        z = 0

        # Grab spherical coordinates of spin vector
        theta = arccos(spin[2])  # theta = arccos(z) since spin is unit vector
        if spin[0] == 0:
            if spin[1] > 0:
                phi = pi/2.0
            else:
                phi = 3*pi/2.0
        else:
            phi = arctan(spin[1]/spin[0]) # phi = arctan(y/x)

        # Rotate random point with phi
        x, y = x*cos(phi) - y*sin(phi), x*sin(phi) + y*cos(phi)

        # Now, rotate random point with theta - rotate around y' = -sin(phi), cos(phi) axis
        u_x = -sin(phi)
        u_y = cos(phi)

        R_matrix = [[cos(theta) + power(u_x,2)*(1-cos(theta)), u_x*u_y*(1-cos(theta)), u_y*sin(theta)],
                    [u_x*u_y*(1-cos(theta)), cos(theta) + power(u_y,2)*(1-cos(theta)), -u_x*sin(theta)],
                    [-u_y*sin(theta), u_x*sin(theta), cos(theta)]]

        final_point = dot(R_matrix, [x, y, z]) + spin
        final_point = final_point/norm(final_point)
        return array(final_point, float) # Return the perturbed spin

    cpdef void simulate(self, int num_sweeps, float T, filename=None):
        cdef int total_accept
        self.sweep_num = 1
        total_accept = 0
        print("Simulating ", num_sweeps, " sweeps at T=", T)
        for i in range(num_sweeps):
            total_accept += self.sweep(T)
        print("Acceptance ratio is: ", total_accept/(num_sweeps*self.lat_size))


    cdef int sweep(self, float T):
        # Perform as many steps as there are lattice sites
        cdef double [:] temp_spin
        cdef double [:] spin
        cdef double [:] delta_spin
        cdef double [:] neighbor_sum
        cdef double delta_a
        cdef int chosen_site_lay, chosen_site_row, chosen_site_col
        cdef int num_accept
        cdef float delta_E
        num_accept = 0

        # Vectorize this loop and move calc_delta_E into this function
        for step in range(2*self.lat_size):
            # Select a new state by randomly choosing a spin to perturb
            # Note: NumPy randint arguments are on a half-open interval
            chosen_site_lay = int(randint(0, 2))
            chosen_site_row = int(randint(0, self.rows))
            chosen_site_col = int(randint(0, self.cols))
            spin = self.lattice[chosen_site_lay][chosen_site_row][chosen_site_col]
            # Calculate the difference in energy between new and old state
            # Using the summation trick of Newman, Barkema (equation 3.10)
            temp_spin = self.perturb(spin)


            #### Moving delta_E calculation into function here ####

            #delta_E = self.calc_delta_E(temp_spin, chosen_site_lay, chosen_site_row, chosen_site_col)

            spin = array(self.lattice[chosen_site_lay][chosen_site_row][chosen_site_col], float)

            delta_spin = array([temp_spin[i] - spin[i] for i in range(3)], float)
            neighbor_sum = self.lattice[chosen_site_lay][(chosen_site_row-1)%self.rows][chosen_site_col] # north neighbor
            neighbor_sum = add(neighbor_sum, self.lattice[chosen_site_lay][(chosen_site_row+1)%self.rows][chosen_site_col]) # south neighbor
            neighbor_sum = add(neighbor_sum, self.lattice[chosen_site_lay][chosen_site_row][(chosen_site_col-1)%self.cols]) # west neighbor
            neighbor_sum = add(neighbor_sum, self.lattice[chosen_site_lay][chosen_site_row][(chosen_site_col+1)%self.cols]) # east neighbor

            if chosen_site_lay == 0:
                delta_a = self.k1*(power(delta_spin[2],2) + 2*delta_spin[2]*spin[2])
            else:
                delta_a = self.k2*(power(delta_spin[2],2) + 2*delta_spin[2]*spin[2])

            delta_E = -self.J_intra*dot(delta_spin, neighbor_sum) + self.J_inter*dot(delta_spin, self.lattice[(chosen_site_lay+1)%2][chosen_site_row][chosen_site_col]) + delta_a - self.B*(delta_spin[2])
            ##### End delta_E calculation ######



            if not ((delta_E > 0) and (rand() >= exp(-(1.0/T)*delta_E))):
                #accept_flag = False
                num_accept += 1
                self.lattice[chosen_site_lay][chosen_site_row][chosen_site_col] = temp_spin
                #self.energy = self.energy + delta_E
                #self.mag = float(self.mag) + 2*float(self.lattice[chosen_site_row][chosen_site_col])

        # After sweep has completed, record energy and magnetization per spin
        #self.energy_vals.append(self.energy)
        #self.mag_vals.append(self.mag/self.lat_size)
        self.sweep_num += 1
        return num_accept
    '''
    cdef float calc_delta_E(self, double [:] temp_spin, int layer, int row, int col):
        # Change in energy given by -J*(delta_spin)*(neighbors)
        cdef double [:] spin
        cdef double [:] delta_spin
        cdef double [:] neighbor_sum
        cdef double delta_a

        spin = array(self.lattice[layer][row][col], float)

        delta_spin = array([temp_spin[i] - spin[i] for i in range(3)], float)

        # Add up all neighbor spins, element wise
        neighbor_sum = self.lattice[layer][(row-1)%self.rows][col] # north neighbor
        neighbor_sum = add(neighbor_sum, self.lattice[layer][(row+1)%self.rows][col]) # south neighbor
        neighbor_sum = add(neighbor_sum, self.lattice[layer][row][(col-1)%self.cols]) # west neighbor
        neighbor_sum = add(neighbor_sum, self.lattice[layer][row][(col+1)%self.cols]) # east neighbor

        if layer == 0:
            delta_a = self.k1*(power(delta_spin[2],2) + 2*delta_spin[2]*spin[2])
        else:
            delta_a = self.k2*(power(delta_spin[2],2) + 2*delta_spin[2]*spin[2])

        return -self.J_intra*dot(delta_spin, neighbor_sum) + self.J_inter*dot(delta_spin, self.lattice[(layer+1)%2][row][col]) + delta_a - self.B*(delta_spin[2])
    '''
    def visualize_lattice(self):
        norm = colors.Normalize(vmin=-1, vmax=1)
        plt.subplot(121)
        plt.imshow([[item[2] for item in row] for row in self.lattice[0]], cmap=plt.get_cmap('Spectral'), norm=norm)
        plt.subplot(122)
        plt.imshow([[item[2] for item in row] for row in self.lattice[1]], cmap=plt.get_cmap('Spectral'), norm=norm)
        plt.show()
    # make d sufficiently large - half j
    # run simulation at different temp - should see spiral phase.
    # plot sz spiral phase cos theta = 1 - theta^2 , sin theta
    # minimize energy and j*theta^2 - d*theta


    # plot m vs h, m1 and m2 v h
    # DM vector - bulk vector pointing from i to j
    # symmetry breaking rij cross z

    cpdef float calc_magnetization(self, int layer):
        ''' Calculate the magnetization of the system, either as a whole or
            for a specific layer.

            :param layer: indicates the layer for which magnetization should be calculated. Defaults
            to None if no layer is selected and all spins are taken into account.
            :type layer: int.
        '''
        cdef float mag, mag_spin
        mag = 0.0
        mag_spin = 0.0
        #print("In calc mag")
        if layer == -1:
            #print("In calc mag -1")

            # Calc magnetization for both layers
            for lay in self.lattice:
                #print("In layer loop, layer = ", lay)
                for row in lay:
                    #print("in row loop, row = ", row)
                    for spin in row:
                        #print("In spin loop, spin = ", spin)
                        mag += spin[2]
            #print("Reassigning mag spin")
            mag_spin = (mag/(self.lat_size*2))
        else:
            #print("In calc mag else")

            # Only magnetization for the first or second layer
            for row in self.lattice[layer]:
                for spin in row:
                    mag += spin[2]
            mag_spin = mag/self.lat_size

        #print("Magnetization is: ", mag)
        print("Magnetization per spin is: ", mag_spin)

        return mag_spin


    def visualize_spins(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make the grid
        for g, layer in enumerate(self.lattice):
            for i, row in enumerate(layer):
                for j, spin in enumerate(row):
                    ax.quiver(i,j,g, spin[0], spin[1], spin[2], length=.4, pivot = 'middle', normalize=True)
        ax.set_zlim(-1,2)
        plt.show()

    def mag_v_temp(self,  init_temp = 0, final_temp = 1, temp_step = .1, eq_time = None, cor_time = None, filename = None):
        ''' Plots the mean magnetization per spin vs temperature for
            lattice. Plot begins with T=0 configuration and plots to final_temp
            in temp_step intervals.
        '''
        mag_vals = []
        mag_vals_1 = []
        mag_vals_2 = []
        temp_vals = []
        curr_temp = init_temp

        while curr_temp > final_temp:
            curr_temp -= temp_step
            self.simulate(eq_time, T = curr_temp)

            self.simulate(cor_time, T = curr_temp)

            avg_mag_per_spin = self.calc_magnetization(-1)  # Mag
            avg_mag_per_spin_1 = self.calc_magnetization(0)  # Mag1
            avg_mag_per_spin_2 = self.calc_magnetization(1)  # Mag2

            mag_vals.append(avg_mag_per_spin)
            mag_vals_1.append(avg_mag_per_spin_1)
            mag_vals_2.append(avg_mag_per_spin_2)

            temp_vals.append(curr_temp)
            print("Current temp: ", curr_temp)

        plt.subplot(311)
        plt.plot(temp_vals, mag_vals, 'go')
        plt.ylabel('$M$')
        plt.xlabel('T')

        plt.subplot(312)
        plt.plot(temp_vals, mag_vals_1, 'go')
        plt.ylabel('$M_{1}$')
        plt.xlabel('T')

        plt.subplot(313)
        plt.plot(temp_vals, mag_vals_2, 'go')
        plt.ylabel('$M_{2}$')
        plt.xlabel('T')
        plt.show()


    cpdef void cool_lattice(self, float T = .5):
        ''' Cool lattice to T = .1 '''
        cdef float curr_temp
        cdef int num_sweeps
        curr_temp = self.init_T
        while curr_temp > T:
            curr_temp -= .045
            print("Cooling to ", curr_temp)

            self.simulate(num_sweeps = 5000, T = curr_temp)

        #with open("h3d_.5_20.txt", 'w') as f:
        #    json.dump(self.lattice, f)

    def record_lattice(self, filename = 'save_lattice.p'):
        l = [[[ list(self.lattice[g][i][j]) for j in range(self.cols)] for i in range(self.rows)] for g in range(2)]
        pickle.dump( l, open( filename, "wb" ) )







def sweep_k():
    cdef float k
    cdef double [:] k_vals
    #cdef float [:] total_mags_per_spin
    #cdef float [:] total_mags_per_spin_1
    #cdef float [:] total_mags_per_spin_2
    start = time.time()
    k_vals = linspace(-5,0,num=100)
    B_SWITCH = .75
    total_mags_per_spin = []
    total_mags_per_spin_1 = []
    total_mags_per_spin_2 = []
    lat = Heisenberg3D(10, 10, k1=-5, k2=-5, J_inter = .5, init_T = 5, B=B_SWITCH)
    #lat.simulate(num_sweeps = 10000, T= .5)
    lat.record_lattice('save_lattice_test.p')
    lat.cool_lattice(.1)
    lat.record_lattice('save_lattice.p')

    for k in k_vals:
        lat.k1 = k
        lat.k2 = k
        print("k = ", k)
        lat.simulate(num_sweeps = 7000, T= .1)
        total_mag_per_spin = lat.calc_magnetization(-1)
        t_mag_per_spin_1 = lat.calc_magnetization(0)
        t_mag_per_spin_2 = lat.calc_magnetization(1)
        total_mags_per_spin.append(total_mag_per_spin)
        total_mags_per_spin_1.append(t_mag_per_spin_1)
        total_mags_per_spin_2.append(t_mag_per_spin_2)

    for k in reversed(k_vals):
        lat.k1 = k
        lat.k2 = k
        print("k = ", k)
        lat.simulate(num_sweeps = 7000, T= .1)
        total_mag_per_spin = lat.calc_magnetization(-1)
        t_mag_per_spin_1 = lat.calc_magnetization(0)
        t_mag_per_spin_2 = lat.calc_magnetization(1)
        total_mags_per_spin.append(total_mag_per_spin)
        total_mags_per_spin_1.append(t_mag_per_spin_1)
        total_mags_per_spin_2.append(t_mag_per_spin_2)

    end = time.time()
    print("Execution time: ", end - start)

    plt.subplot(311)   # Magnetization per spin, both lattices
    plt.plot(list(k_vals), total_mags_per_spin[0:len(list(k_vals))], 'r')
    plt.plot(list(reversed(k_vals)), total_mags_per_spin[len(list(k_vals)):], 'b')

    plt.ylabel("M")
    plt.xlabel("k")

    #plt.subplot(323)   # Total magnetization, first lattice
    #plt.plot(k_vals, t_mags_1)

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(list(k_vals), total_mags_per_spin_1[0:len(list(k_vals))], 'r')
    plt.plot(list(reversed(k_vals)), total_mags_per_spin_1[len(list(k_vals)):], 'b')
    plt.ylabel("$M_{1}$")
    plt.xlabel("k")


    #plt.subplot(325)   # Total magnetization, second lattice
    #plt.plot(k_vals, t_mags_2)

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(list(k_vals), total_mags_per_spin_2[0:len(list(k_vals))], 'r')
    plt.plot(list(reversed(k_vals)), total_mags_per_spin_2[len(list(k_vals)):], 'b')
    plt.ylabel("$M_{2}$")
    plt.xlabel("k")
    plt.show()






def plot_M_v_B():
    logging.basicConfig(filename="bilayer_h.log", filemode='w',level=logging.INFO, format='%(message)s')
    print("50,000 sweeps per measure.")
    total_mags_per_spin = []
    total_mags_per_spin_1 = []
    total_mags_per_spin_2 = []
    B_vals = linspace(-5,5,num=100)
    lat = Heisenberg3D(20, 20, k1=-1, k2=-1, J_inter = 1, init_T = 5, B=B_vals[0])
    #lat.simulate(num_sweeps = 10000, T= .5)
    lat.cool_lattice(.5)

    #print("B_vals: ", B_vals)
    #print("reversed: ", list(reversed(B_vals)))
    #print("added: ", list(B_vals)+list(reversed(B_vals)))
    logging.info('**************************')
    logging.info('Creating a 3D lattice with ' + str(lat.lat_size) + ' sites per layer.')
    logging.info('B: ' + str(lat.B))
    logging.info('M: ' + str())
    for B in B_vals:
        # Make k negative   --- flipped k value, bad effects
        # Play with temp    - Cooling lattice slowly.
        # Mag v temp        - plotted mag v temp. Results as expected
        # strong interlayer coupling     - makes results more pronounced, but not correct
        # sweeping speed     - Reduced to 1000 sweeps per measurement
        # Should have opposite magnetization on layers - investigate

        print("B: ", B)
        lat.B = B
        lat.simulate(num_sweeps = 5000, T= .5)
        total_mag, total_mag_per_spin = lat.calc_magnetization(-1)
        t_mag_1, t_mag_per_spin_1 = lat.calc_magnetization(0)
        t_mag_2, t_mag_per_spin_2 = lat.calc_magnetization(1)
        #total_mags.append(total_mag)
        total_mags_per_spin.append(total_mag_per_spin)
        #t_mags_1.append(t_mag_1)
        #t_mags_2.append(t_mag_2)
        total_mags_per_spin_1.append(t_mag_per_spin_1)
        total_mags_per_spin_2.append(t_mag_per_spin_2)
        logging.info('B: ' + str(B))
        logging.info('M: ' + str(total_mag_per_spin))

    for B in reversed(B_vals):
        # Make k negative
        print("B: ", B)
        lat.B = B
        lat.simulate(num_sweeps = 5000, T= .5)
        total_mag, total_mag_per_spin = lat.calc_magnetization(-1)
        t_mag_1, t_mag_per_spin_1 = lat.calc_magnetization(0)
        t_mag_2, t_mag_per_spin_2 = lat.calc_magnetization(1)
        #total_mags.append(total_mag)
        total_mags_per_spin.append(total_mag_per_spin)
        #t_mags_1.append(t_mag_1)
        #t_mags_2.append(t_mag_2)
        total_mags_per_spin_1.append(t_mag_per_spin_1)
        total_mags_per_spin_2.append(t_mag_per_spin_2)
        logging.info('B: ' + str(B))
        logging.info('M: ' + str(total_mag_per_spin))
    #plt.subplot(321)   # Total magnetization, both lattices
    #plt.plot(k_vals, total_mags)

    plt.subplot(311)   # Magnetization per spin, both lattices
    plt.plot(list(B_vals), total_mags_per_spin[0:len(list(B_vals))], 'r')
    plt.plot(list(reversed(B_vals)), total_mags_per_spin[len(list(B_vals)):], 'b')

    plt.ylabel("M")
    plt.xlabel("B")

    #plt.subplot(323)   # Total magnetization, first lattice
    #plt.plot(k_vals, t_mags_1)

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(list(B_vals), total_mags_per_spin_1[0:len(list(B_vals))], 'r')
    plt.plot(list(reversed(B_vals)), total_mags_per_spin_1[len(list(B_vals)):], 'b')
    plt.ylabel("$M_{1}$")
    plt.xlabel("B")


    #plt.subplot(325)   # Total magnetization, second lattice
    #plt.plot(k_vals, t_mags_2)

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(list(B_vals), total_mags_per_spin_2[0:len(list(B_vals))], 'r')
    plt.plot(list(reversed(B_vals)), total_mags_per_spin_2[len(list(B_vals)):], 'b')
    plt.ylabel("$M_{2}$")
    plt.xlabel("B")
    plt.show()

def plot_M_v_k(B = 0):
    total_mags = []
    total_mags_per_spin = []
    t_mags_1 = []
    total_mags_per_spin_1 = []
    t_mags_2 = []
    total_mags_per_spin_2 = []
    k_vals = linspace(-1,0,num=100)
    for k in k_vals:
        lat = Heisenberg3D(10,10, k1=k, k2=k,init_T = 2, B = B)
        lat.simulate(num_sweeps = 5000, T= .001)
        total_mag, total_mag_per_spin = lat.calc_magnetization(-1)
        t_mag_1, t_mag_per_spin_1 = lat.calc_magnetization(0)
        t_mag_2, t_mag_per_spin_2 = lat.calc_magnetization(1)
        total_mags.append(total_mag)
        total_mags_per_spin.append(total_mag_per_spin)
        t_mags_1.append(t_mag_1)
        t_mags_2.append(t_mag_2)
        total_mags_per_spin_1.append(t_mag_per_spin_1)
        total_mags_per_spin_2.append(t_mag_per_spin_2)
    #plt.subplot(321)   # Total magnetization, both lattices
    #plt.plot(k_vals, total_mags)

    plt.subplot(311)   # Magnetization per spin, both lattices
    plt.plot(k_vals, total_mags_per_spin)
    plt.ylabel("Magnetization per spin - entire")
    plt.xlabel("k strength")

    #plt.subplot(323)   # Total magnetization, first lattice
    #plt.plot(k_vals, t_mags_1)

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(k_vals, total_mags_per_spin_1)
    plt.ylabel("Magnetization per spin - lattice 1")
    plt.xlabel("k strength")


    #plt.subplot(325)   # Total magnetization, second lattice
    #plt.plot(k_vals, t_mags_2)

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(k_vals, total_mags_per_spin_2)
    plt.ylabel("Magnetization per spin - lattice 2")
    plt.xlabel("k strength")



    plt.show()







if __name__ == "__main__":
    start = time.time()
    #lat = Heisenberg3D(10,10, init_T = 5)
    #lat.simulate(num_sweeps = 10000, T=.001)
    #lat.calc_magnetization()
    #lat.calc_magnetization(0)
    #lat.calc_magnetization(1)

    #lat.mag_v_temp(init_temp = 5, final_temp=0.01, temp_step=.05, eq_time=7000, cor_time=1000)

    #plot_M_v_B()
    sweep_k()
    #plot_M_v_k()
    end = time.time()

    print("Execution time (cython): ", end - start)
