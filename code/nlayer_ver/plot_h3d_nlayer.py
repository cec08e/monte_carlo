from __future__ import print_function
import ctypes
import pickle

SIM_NUM = 10
NUM_L = 3

PDOUBLE = ctypes.POINTER(ctypes.c_double)
PPDOUBLE = ctypes.POINTER(PDOUBLE)

_h3d = ctypes.CDLL('./heisenberg3d_nlayer.so')
_h3d.initialize_lattice.argtypes = ()
_h3d.M_v_B.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.M_v_B.restype = ctypes.c_int
#_h3d.M_v_T.argtypes = (ctypes.POINTER(PDOUBLE), ctypes.c_double, ctypes.c_double, ctypes.c_double)
#_h3d.M_v_T.restype = ctypes.c_int
#_h3d.M_v_K.argtypes = (ctypes.POINTER(PDOUBLE),)
#_h3d.M_v_K.restype = ctypes.c_int
#_h3d.M_v_J.argtypes = (ctypes.POINTER(PDOUBLE),)
#_h3d.M_v_J.restype = ctypes.c_int

def plot_M_v_B(max_samples = 10000):
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
    sample_arr = ctypes.c_double*(NUM_L + 2) # B, M, M1, M2, M3, M4, ... array for specific B value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()

    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_B(results) # pass results list by ref to be populated by C function
    ''''
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
    '''
    data = [[] for i in range((NUM_L + 2))]
    for i in range(num_samples):
        for j in range(NUM_L + 2):
            data[j].append(results[i][j])
    '''
        B_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])
        M3_vals.append(results[i][4])
        M4_vals.append(results[i][5])
    '''
    '''
    plt.subplot(511)   # Magnetization per spin, all lattices
    plt.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)

    plt.ylabel("M")
    #plt.xlabel("B")
    plt.xticks([])

    plt.subplot(512)   # Magnetization per spin, first lattice
    plt.plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{1}$")
    plt.xticks([])
    #plt.xlabel("B")

    plt.subplot(513)   # Magnetization per spin, second lattice
    plt.plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{2}$")
    plt.xticks([])

    #plt.xlabel("B")

    plt.subplot(514)   # Magnetization per spin, third lattice
    plt.plot(B_vals[0:int(num_samples/2)], M3_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M3_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{3}$")
    #plt.xlabel("B")
    plt.xticks([])


    plt.subplot(515)   # Magnetization per spin, fourth lattice
    plt.plot(B_vals[0:int(num_samples/2)], M4_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M4_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{4}$")
    plt.xlabel("B (J_inter = {.05,.05,.05,0}, K = {.05,.05,.05,.05}, J_intra = {1.0,1.0,1.0,1.0})")
    #plt.xlabel("$H$ ($J_{FM} = 1.0$, $J_{AFM} = .05$, $K_{1} = .06$, $K_{234} = .05$)")
    plt.savefig("sim_results/sim_"+str(SIM_NUM)+".png")
    '''
    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)

    #plt.show()

if __name__ == "__main__":
    plot_M_v_B()
