from __future__ import print_function
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors
import ctypes

_h3d = ctypes.CDLL('heisenberg3d.so')
_h3d.initialize_lattice.argtypes = ()
_h3d.M_v_B.argtypes = ()
_h3d.mag_v_temp.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int)

def mag_v_temp(init_temp = 0, final_temp = 1, temp_step = .1, eq_time = None, cor_time = None):
    ''' Plots the mean magnetization per spin vs temperature for
        lattice. Plot begins with T=0 configuration and plots to final_temp
        in temp_step intervals.
    '''
    mag_vals = []
    mag_vals_1 = []
    mag_vals_2 = []
    temp_vals = []
    global _h3d
    print("Initializing lattice...")
    _h3d.initialize_lattice()
    results = _h3d.mag_v_temp(5.0, .01, .05, 10000, 5000)

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

def plot_M_v_B():
    global _h3d
    print("Initializing lattice...")
    _h3d.initialize_lattice()
    results = _h3d.M_v_B()
    print(results)
    print(type(results))
    '''
    plt.subplot(311)   # Magnetization per spin, both lattices
    plt.plot(list(B_vals), total_mags_per_spin[0:len(list(B_vals))], 'r')
    plt.plot(list(reversed(B_vals)), total_mags_per_spin[len(list(B_vals)):], 'b')

    plt.ylabel("M")
    plt.xlabel("B")

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(list(B_vals), total_mags_per_spin_1[0:len(list(B_vals))], 'r')
    plt.plot(list(reversed(B_vals)), total_mags_per_spin_1[len(list(B_vals)):], 'b')
    plt.ylabel("$M_{1}$")
    plt.xlabel("B")

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(list(B_vals), total_mags_per_spin_2[0:len(list(B_vals))], 'r')
    plt.plot(list(reversed(B_vals)), total_mags_per_spin_2[len(list(B_vals)):], 'b')
    plt.ylabel("$M_{2}$")
    plt.xlabel("B")
    plt.show()
    '''


if __name__ == "__main__":
    #mag_v_temp(5, .1, 10000, 5000)
    plot_M_v_B()
