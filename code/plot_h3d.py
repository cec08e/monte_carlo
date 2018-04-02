from __future__ import print_function
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors
import ctypes


PDOUBLE = ctypes.POINTER(ctypes.c_double)
PPDOUBLE = ctypes.POINTER(PDOUBLE)


_h3d = ctypes.CDLL('heisenberg3d.so')
_h3d.initialize_lattice.argtypes = ()
_h3d.M_v_B.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.M_v_B.restype = ctypes.c_int
_h3d.mag_v_temp.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int)
_h3d.M_v_K.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.M_v_K.restype = ctypes.c_int

'''
class Spin(ctypes.Structure):
    _fields_ = [('x', ctypes.c_double), ('y', ctypes.c_double), ('z', ctypes.c_double)]
'''


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

def plot_M_v_B(max_samples = 500):
    global _h3d

    # Create results list to pass by reference
    sample_arr = ctypes.c_double*4 # B, M, M1, M2 array for specific B value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()


    print("Initializing lattice...")
    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_B(results) # pass results by ref

    B_vals = []
    M_vals = []
    M1_vals = []
    M2_vals = []
    for i in range(num_samples):
        B_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])


    plt.subplot(311)   # Magnetization per spin, both lattices
    plt.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')


    plt.ylabel("M")
    plt.xlabel("B")

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')

    plt.ylabel("$M_{1}$")
    plt.xlabel("B")

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')

    plt.ylabel("$M_{2}$")
    plt.xlabel("B")

    plt.show()

def plot_M_v_K(max_samples = 500):
    global _h3d

    # Create results list to pass by reference
    sample_arr = ctypes.c_double*4 # B, M, M1, M2 array for specific B value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()


    print("Initializing lattice...")
    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_K(results) # pass results by ref

    K_vals = []
    M_vals = []
    M1_vals = []
    M2_vals = []
    for i in range(num_samples):
        K_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])


    plt.subplot(311)   # Magnetization per spin, both lattices
    plt.plot(K_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    plt.plot(K_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')


    plt.ylabel("M")
    plt.xlabel("K")

    plt.subplot(312)   # Magnetization per spin, first lattice
    plt.plot(K_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    plt.plot(K_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')

    plt.ylabel("$M_{1}$")
    plt.xlabel("K")

    plt.subplot(313)   # Magnetization per spin, second lattice
    plt.plot(K_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    plt.plot(K_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')

    plt.ylabel("$M_{2}$")
    plt.xlabel("K")

    plt.show()

if __name__ == "__main__":
    #mag_v_temp(5, .1, 10000, 5000)
    plot_M_v_B()
    #plot_M_v_K()
