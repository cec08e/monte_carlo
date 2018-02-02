from __future__ import print_function
from numpy.random import uniform
from numpy.random import random, randint, rand
from numpy import power,sqrt,pi,sin,cos,arccos,arctan
from numpy import dot, cross, add, exp
from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors



class Heisenberg3D(object):
    ''' Bilayer implementation '''
    def __init__(self, rows, cols, init_T = 0, B = 0, J_intra = 1, J_inter = .1, k1 = 0, k2 = 0, D = (1,0,0)):
        self.rows = rows
        self.cols = cols
        self.J_intra = J_intra     # Intra-layer interaction
        self.J_inter = J_inter     # Inter-layer interaction
        self.D = D
        self.k1 = k1
        self.k2 = k2
        self.lat_size = rows*cols
        self.initialize(init_T, B)
        self.sweep_num=0

    def initialize(self, init_T, B):
        '''
        Initializes lattice with spin configuration
        determined by temperature.
        '''
        # Generate random spin configuration
        # Two layers
        self.lattice = [[[ self.gen_random_spin() for j in range(self.cols)] for i in range(self.rows)] for g in range(2)]


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
        - Generate point x,y in circle with radius R
          see: http://xdpixel.com/random-points-in-a-circle/
        - Place point on tangent plane defined by current spin vector
            - Shift x,y origin by adding current spin vector.
            - Rotate plane by theta of current spin vector.

        '''
        R = 1
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(*spin, c='r', marker='^')
        #print("Spin vector is ", spin)
        #const = (power(spin[0],2) + power(spin[1],2) + power(spin[2],2))
        #print("Tangent plane defined by const = ", const)
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
        #print(final_point)
        #print("Equation result: ", (final_point[0]*spin[0] + final_point[1]*spin[1] + final_point[2]*spin[2]))
        #ax.scatter(*final_point, c='b', marker='o')
        #ax.set_xlim(-1,1)
        #ax.set_ylim(-1,1)
        #ax.set_zlim(-1,1)
        #plt.show()
        return tuple(final_point) # Return the perturbed spin

    def simulate(self, num_sweeps, T, filename=None):
        self.sweep_num = 1
        for i in range(num_sweeps):
            self.sweep(T)
            if i+1 in [1,10,100,1000,10000]:
                self.visualize_lattice()
                self.visualize_spins()


    def sweep(self, T):
        # Perform as many steps as there are lattice sites
        for step in range(self.lat_size):
            # Select a new state by randomly choosing a spin to perturb
            # Note: NumPy randint arguments are on a half-open interval
            chosen_site_lay = randint(0, 2)
            chosen_site_row = randint(0, self.rows)
            chosen_site_col = randint(0, self.cols)
            spin = self.lattice[chosen_site_lay][chosen_site_row][chosen_site_col]
            # Calculate the difference in energy between new and old state
            # Using the summation trick of Newman, Barkema (equation 3.10)
            temp_spin = self.perturb(spin)
            delta_E = self.calc_delta_E(temp_spin, chosen_site_lay, chosen_site_row, chosen_site_col)
            #print("Delta E: ", delta_E)

            if not ((delta_E > 0) and (rand() >= exp(-(1.0/T)*delta_E))):   # FIX THISSSSS
                #accept_flag = False
                self.lattice[chosen_site_lay][chosen_site_row][chosen_site_col] = temp_spin
                #self.energy = self.energy + delta_E
                #self.mag = float(self.mag) + 2*float(self.lattice[chosen_site_row][chosen_site_col])

        # After sweep has completed, record energy and magnetization per spin
        #self.energy_vals.append(self.energy)
        #self.mag_vals.append(self.mag/self.lat_size)
        self.sweep_num += 1

    def calc_delta_E(self, temp_spin, layer, row, col):
        # Change in energy given by -J*(delta_spin)*(neighbors)
        spin = self.lattice[layer][row][col]
        delta_spin = [temp_spin[i] - spin[i] for i in range(3)]

        # Add up all neighbor spins, element wise
        neighbor_sum = self.lattice[layer][(row-1)%self.rows][col] # north neighbor
        neighbor_sum = add(neighbor_sum, self.lattice[layer][(row+1)%self.rows][col]) # south neighbor
        neighbor_sum = add(neighbor_sum, self.lattice[layer][row][(col-1)%self.cols]) # west neighbor
        neighbor_sum = add(neighbor_sum, self.lattice[layer][row][(col+1)%self.cols]) # east neighbor

        #cross_term = cross(delta_spin, neighbor_sum)
        if layer == 0:
            delta_a = self.k1*power(delta_spin[2],2)
        else:
            delta_a = self.k2*power(delta_spin[2],2)

        return -self.J_intra*dot(delta_spin, neighbor_sum) + self.J_inter*dot(delta_spin, self.lattice[(layer+1)%2][row][col]) + delta_a

    def visualize_lattice(self):
        norm = colors.Normalize(vmin=-1, vmax=1)
        #cmap = colors.ListedColormap(['white','black'])
        #plt.imshow([[item[2] for item in row] for row in self.lattice], cmap=cmap, norm=norm)
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


if __name__ == "__main__":
    lat = Heisenberg3D(20,20, init_T = 2)
    lat.simulate(num_sweeps = 10000, T=.001)
