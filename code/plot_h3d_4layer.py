from __future__ import print_function
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors
import ctypes

"""plot_h3d_4layer.py

This module provides plotting functionality and controls the calling of
MC simulation functions in the associated C program. Currently, supported
plotting routines include Magnetization vs External Field (Anisotropy and
interlayer coupling fixed), Magnetization vs Anisotropy (external field
and interlaying coupling fixed), and Magnetization vs Temperature (all
else held fixed).

Requirements:
    - Python 3
    - matplotlib
    - ctypes (standard)

Example:
        $ python3 plot_h3d_4layer.py
-- Currently: Make sure to update the desired call in the main conditional at
the bottom of the module.
-- Make sure to update labels on graphs when changing the simulation parameters.

Attributes:
    This module makes use of the ctypes library to communicate with the MC
    code that is written in C. To that end, there are a number of top-level
    module attributes such as PDOUBLE, _h3d, etc. However, a user of this
    module should not need to use and/or modify these attributes.

Todo:
    * Add command line flags to call module functions.

"""
SIM_NUM = 20

PDOUBLE = ctypes.POINTER(ctypes.c_double)
PPDOUBLE = ctypes.POINTER(PDOUBLE)

_h3d = ctypes.CDLL('heisenberg3d_4layer.so')
_h3d.initialize_lattice.argtypes = ()
_h3d.M_v_B.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.M_v_B.restype = ctypes.c_int
_h3d.M_v_T.argtypes = (ctypes.POINTER(PDOUBLE), ctypes.c_double, ctypes.c_double, ctypes.c_double)
_h3d.M_v_T.restype = ctypes.c_int
_h3d.M_v_K.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.M_v_K.restype = ctypes.c_int

def plot_M_v_T(init_temp = 0, final_temp = 1, temp_step = .1):
    ''' Plots the mean magnetization per spin vs temperature for
        lattice. Plot begins with T=0 configuration and plots to final_temp
        in temp_step intervals.
    '''
    max_samples = 5000
    M_vals = []
    M1_vals = []
    M2_vals = []
    M3_vals = []
    M4_vals = []
    temp_vals = []
    global _h3d

    sample_arr = ctypes.c_double*6 # T, M, M1, M2, M3, M4 array for specific T value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()

    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_T(results, 5.0, .01, .05)

    for i in range(num_samples):
        temp_vals.append(results[i][0])
        M_vals.append(results[i][4])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])
        M3_vals.append(results[i][4])
        M4_vals.append(results[i][5])

    plt.subplot(311)
    plt.plot(temp_vals, M_vals, 'g')
    plt.ylabel('$M1 - M2$')
    plt.xlabel('T')

    plt.subplot(312)
    plt.plot(temp_vals, M1_vals, 'g')
    plt.ylabel('$M_{1}$')
    plt.xlabel('T')

    plt.subplot(313)
    plt.plot(temp_vals, M2_vals, 'g')
    plt.ylabel('$M_{2}$')
    plt.xlabel('T')

    plt.subplot(312)
    plt.plot(temp_vals, M3_vals, 'g')
    plt.ylabel('$M_{3}$')
    plt.xlabel('T')

    plt.subplot(312)
    plt.plot(temp_vals, M4_vals, 'g')
    plt.ylabel('$M_{4}$')

    plt.xlabel('T (J_inter = {.1, .1, .1, 0} , K = {.1, .1, .1, .1})')
    plt.show()

def plot_M_v_B(max_samples = 5000):
    """ Plots the mean magnetization per spin vs. external field strength.
        Data is acquired by calling the M_v_B() function in the C implementation.

        Plots will include separate sub-figures for total magnetization M, and
        individual layer magnetizations (identified as M1, M2, M3, M4).

        The number of data points returned by the C function is available as
        num_samples. The data will populate the results list, where results[i][0]
        is the external field strength of the data point and results[i][1] is the
        total magnetization. The data results[i][2:5] represent the magnetization
        values of the individual layers.
    """
    global _h3d

    # Create results list to pass by reference
    sample_arr = ctypes.c_double*6 # B, M, M1, M2, M3, M4 array for specific B value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()

    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_B(results) # pass results list by ref to be populated by C function

    B_vals = []
    M_vals = []
    M1_vals = []
    M2_vals = []
    M3_vals = []
    M4_vals = []
    for i in range(num_samples):
        B_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])
        M3_vals.append(results[i][4])
        M4_vals.append(results[i][5])


    plt.subplot(511)   # Magnetization per spin, all lattices
    plt.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)

    plt.ylabel("M")
    plt.xlabel("B")

    plt.subplot(512)   # Magnetization per spin, first lattice
    plt.plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{1}$")
    plt.xlabel("B")

    plt.subplot(513)   # Magnetization per spin, second lattice
    plt.plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{2}$")
    plt.xlabel("B")

    plt.subplot(514)   # Magnetization per spin, third lattice
    plt.plot(B_vals[0:int(num_samples/2)], M3_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M3_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{3}$")
    plt.xlabel("B")

    plt.subplot(515)   # Magnetization per spin, fourth lattice
    plt.plot(B_vals[0:int(num_samples/2)], M4_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M4_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{4}$")
    plt.xlabel("B (J_inter = {.05,.05,.05,0}, K = {.06,.05,.05,.05}, J_intra = 1)")

    plt.savefig("sim_results/sim_"+str(SIM_NUM)+".png")

    plt.show()




def plot_M_v_K(max_samples = 5000):
    """ Plots the mean magnetization per spin vs. anisotropic strength.
        Data is acquired by calling the M_v_K() function in the C implementation.

        Plots will include separate sub-figures for total magnetization M, and
        individual layer magnetizations (identified as M1, M2, M3, M4).

        The number of data points returned by the C function is available as
        num_samples. The data will populate the results list, where results[i][0]
        is the external field strength of the data point and results[i][1] is the
        total magnetization. The data results[i][2:5] represent the magnetization
        values of the individual layers.
    """
    global _h3d

    # Create results list to pass by reference
    sample_arr = ctypes.c_double*6 # B, M, M1, M2, M3, M4 array for specific K value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()

    print("Initializing lattice...")
    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_K(results) # pass results list by ref to be populated by C function

    K_vals = []
    M_vals = []
    M1_vals = []
    M2_vals = []
    M3_vals = []
    M4_vals = []
    for i in range(num_samples):
        K_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])
        M3_vals.append(results[i][4])
        M4_vals.append(results[i][5])


    plt.subplot(311)   # Magnetization per spin, all lattices
    plt.plot(K_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    plt.plot(K_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("M")
    plt.xlabel("K")

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(K_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    plt.plot(K_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{1}$")
    plt.xlabel("K")

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(K_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    plt.plot(K_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{2}$")
    plt.xlabel("K")

    plt.subplot(514)   # Magnetization per spin, third lattice
    plt.plot(B_vals[0:int(num_samples/2)], M3_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M3_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{3}$")
    plt.xlabel("K")

    plt.subplot(515)   # Magnetization per spin, fourth lattice
    plt.plot(B_vals[0:int(num_samples/2)], M4_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M4_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{4}$")
    plt.xlabel("K (J_inter = {.05,.05,.05,0}, B = .1)")

    plt.show()

if __name__ == "__main__":
    #mag_v_temp(5, .1, 10000, 5000)
    #plot_M_v_T(5, .1, .1)
    plot_M_v_B()
    #plot_M_v_K()
