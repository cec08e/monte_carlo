import pickle
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import cycle
#import matplotlib.gridspec as gridspec

orange = "#FF8000"
olive = "#008000"

paper_data = [[2.5, -0.168000],
                [2.4, -0.053333],
                [2.3, -0.015667],
                [2.2, -0.066600],
                [2.1,  -0.163333],
                [2.0, -0.025800],
                [1.9, -0.171000 ],
                [1.8, -0.291667],
                [1.7, -0.309600],
                [1.6, -0.351667],
                [1.5, -0.244000],
                [1.4,  -0.488000 ],
                [1.3, -0.559000],
                [1.2, -0.617000 ],
                [1.1, -1.085667 ],
                [1.0, -1.251600],
                [0.9, -1.283800],
                [0.8, -1.244200],
                [0.7, -.937400],
                [0.6, -0.409000],
                [0.5, -0.126000 ],
                [0.4, -0.014667],
                [0.3, -0.000667 ],
                [0.2, 0.0]]

paper_data2 = [[2.5, -0.036602],
                [2.4, -0.038283],
                [2.3, -0.049154],
                [2.2, -0.059861 ],
                [2.1,  -0.088860],
                [2.0, -0.097779],
                [1.9, -0.116230 ],
                [1.8, -0.155025],
                [1.7, -0.202004],
                [1.6, -0.255624],
                [1.5, -0.323930],
                [1.4, -0.413181],
                [1.3, -0.556935],
                [1.2, -0.750788],
                [1.1, -0.980471],
                [1.0, -1.224319],
                [0.9, -1.355731],
                [0.8, -1.248796],
                [0.7, -0.868108],
                [0.6, -0.409000],
                [0.5, -0.126000],
                [0.4, -0.014667],
                [0.3, -0.000667],
                [0.2, 0.0]]

paper_data_AF = [[1.800000, -0.011119],[1.600000, -0.019056],[1.400000, -0.034564], [0.800000, -0.028941],
                  [2.200000, -0.006753], [ 0.700000, -0.057552], [1.200000, -0.034000], [ 0.500000, -0.140218],
                  [2.00000,  -0.017915], [1.000000, -0.030311], [.4, -0.253787], [.3, -0.372230]]

paper_data_AF_aux = [ [0,0.006728, 'b'],  [.05, -0.008703, 'r'], [.15, -0.032307, 'c'], [.2,  -0.028941, 'g'], [.25, -0.050739, 'm'], [.3, -0.064178, 'y'], [.35, -0.068622, 'k'], [.5,-0.092457, orange ]]

Q_v_B_F = [[0.1, -0.017641], [0.15,	-0.026679], [0.2,	-0.03724], [0.25,	-0.046759], [0.3,	-0.046819]]

Q_v_B_AF = [[0.1, -0.007179],
[0.125,	0.003636],
[0.15, -0.003991],
[0.175, -0.003750],
[0.2, -0.007241],
[0.225, -0.006181],
[0.25, -0.004236],
[0.275, -0.004767],
[0.3, -0.010002]]

Q_v_B_AF_2 = [[0.1,	-0.007519],
[0.15, -0.006679],
[0.2, -0.017741],
[0.25, -0.02521],
[0.3, -0.024545]]

Q_v_B_AF_3 = [[0.0025, -0.54109],
[0.005, -0.735952],
[0.0075, -0.884773],
[0.01, -0.994004],
[0.0125, -1.065265],
[0.015,	-1.125666],
[0.0175, -1.184708],
[0.02, -1.216218],
[0.0225, -1.251148],
[0.025, -1.269478],
[0.0275,-1.279974],
[0.03, -1.300069],
[0.0325,-1.313777],
[0.035,	-1.307728],
[0.0375,-1.314436],
[0.04, -1.310388],
[0.0425,-1.312315],
[0.045,	-1.304372],
[0.0475,-1.299307],
[0.05,	-1.298029]]

Q_v_B_AF_4 = [
[0.01, -0.001059],
[0.015, -0.008161],
[0.02, -0.024413],
[0.025, -0.049514],
[0.03, -0.078246],
[0.035,	-0.112615],
[0.04, -0.146523],
[0.045, -0.177181],
[0.05, -0.208513],
[0.055, -0.237658],
[0.06, -0.260273],
[0.065, -0.283633],
[0.07, -0.309025],
[0.075, -0.320559],
[0.08, -0.344905],
[0.085, -0.503599],
[0.09, -0.376266],
[0.095, -0.386985],
[0.1, -0.398437],
[0.15, -0.48154],
[0.2, -0.529638],
[0.25, -0.523725],
[0.3, -0.495835],
[0.35, -0.417266]
]

Q_v_B_AF_5 = [
[0.01,0.000001],
[0.03,-0.003717],
[0.05,-0.027315],
[0.07,-0.053695],
[0.09,-0.084335],
[0.11,-0.101177],
[0.13,-0.117126],
[0.15,-0.13034],
[0.17,-0.142905],
[0.19,-0.146039],
[0.21,-0.166408],
[0.23,-0.173344],
[0.25,-0.192012]
]

def visualize_lattice(sim_num):
    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

    data = pickle.load(open(filename, 'rb'))
    print(data)
    plt.imshow([[item[2] for item in row] for row in data[0]], cmap=plt.get_cmap('Spectral'))
    plt.show()

def visualize_spins(sim_num):

    filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

    data = pickle.load(open(filename, 'rb'))
    print(data)
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make the grid
    for i, row in enumerate(data[0]):
        for j, spin in enumerate(row):
            print("Spin: ", spin)
            ax.quiver(i*2, j*2, 0, spin[0], spin[1], spin[2], length=1, pivot = 'tail', normalize=True)
    ax.set_zlim(-30,30)
    plt.show()


def plot_auto_data(data_list):
    col = 'g'
    lab = "0 OR"
    i = 0
    for sim in data_list:
        filename = "sim_results/sim_" + str(sim) + ".pickle"
        data = pickle.load(open(filename, 'rb'))
        print(data)
        plt.plot(data[0], data[1], color=col, label =lab)
        if i == 0:
            col = 'b'
            lab = "1 OR"
        if i == 1:
            col = 'r'
            lab = '2 OR'
        i = 1
        plt.xlim(-10,200)
    plt.xlabel("Sweeps")
    plt.ylabel("$\chi$")
    plt.legend()
    plt.show()


def plot_paper_data():
    #plt.plot([x[0] for x in paper_data_AF], [x[1]*.976 for x in paper_data_AF], 'g^')
    plt.ylabel("Q/1000 spins")
    plt.xlabel("T")
    #plt.ylim(-1.75,.1)
    #plt.xlim(0,3.0)
    #plt.text(1.5, -.8, "J = -1.00, D = .3, B = .2, L = 32", withdash=False)
    #plt.text(1.0, -.1, "J = -1.00, D = .3, T = .8, L = 32", withdash=False)

    '''
    for x in Q_v_B_F:
        #print(x[0], x[1], x[2])
        plt.plot(x[0], x[1], '^r')

    for x in Q_v_B_AF:
        #print(x[0], x[1], x[2])
        plt.plot(x[0], x[1], 'og')

    for x in Q_v_B_AF_2:
        #print(x[0], x[1], x[2])
        plt.plot(x[0], x[1], 'ob')
    '''

    #F_line, = plt.plot([x[0] for x in Q_v_B_F], [x[1]*(1000/400.0) for x in Q_v_B_F], '^r', label="$J=1$, $D=.3$")
    #AF_line_1, = plt.plot([x[0] for x in Q_v_B_AF], [x[1]*(1000/400.0) for x in Q_v_B_AF], 'ob', label="$J=-1$, $D=.3$")
    #AF_line_2, = plt.plot([x[0] for x in Q_v_B_AF_2], [x[1]*(1000/400.0) for x in Q_v_B_AF_2], 'og', label="$J=-1$, $D=.5$")
    AF_line_5, = plt.plot([x[0] for x in Q_v_B_AF_5], [x[1]*(1000/400.0) for x in Q_v_B_AF_5], 'og', label="$D=.3$, $B=1.45$")

    plt.legend(handles=[AF_line_5], title="N = 20")


    plt.show()

def plot_k_data(sim_num):
    filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    data = pickle.load(open(filename, 'rb'))
    print(data)
    num_samples = len(data[0])
    B_vals = data[0]
    M_vals = data[1]
    M1_vals = data[2]
    M2_vals = data[3]
    M3_vals = data[4]
    M4_vals = data[5]

    fig, ax = plt.subplots(5, 1, figsize=(10,10), sharex=True)   # Magnetization per spin, all lattices
    ax[0].plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    ax[0].plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')
    ax[0].set_ylim(-1.0, 1.0)
    ax[0].set_xlim(-1.05,0.025)

    ax[0].set_ylabel("M",fontsize=13)
    #plt.xlabel("B")
    #plt.xticks([])

    #plt.subplot(512)   # Magnetization per spin, first lattice
    ax[1].plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    ax[1].plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')
    ax[1].set_ylim(-1.0, 1.0)


    ax[1].set_ylabel("$M_{1}$",fontsize=13)
    #ax[1].text(.3, -.5, "$K = .2$\n$J_{AFM} = .3$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    #plt.xticks([])
    #plt.xlabel("B")

    #plt.subplot(513)   # Magnetization per spin, second lattice
    ax[2].plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    ax[2].plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    ax[2].set_ylim(-1.0, 1.0)


    ax[2].set_ylabel("$M_{2}$",fontsize=13)
    #ax[2].text(.3, -.5, "$K = .05$\n$J_{AFM} = .05$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))

    #plt.xticks([])

    #plt.xlabel("B")

    #plt.subplot(514)   # Magnetization per spin, third lattice
    ax[3].plot(B_vals[0:int(num_samples/2)], M3_vals[0:int(num_samples/2)], 'r')
    ax[3].plot(B_vals[int(num_samples/2):], M3_vals[int(num_samples/2):], 'b')
    ax[3].set_ylim(-1.0, 1.0)


    ax[3].set_ylabel("$M_{3}$",fontsize=13)
    ax[3].text(.3, -.5, "$K = .05$\n$J_{AFM} = .05$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))

    #plt.xlabel("B")
    #plt.xticks([])


    #plt.subplot(515)   # Magnetization per spin, fourth lattice
    ax[4].plot(B_vals[0:int(num_samples/2)], M4_vals[0:int(num_samples/2)], 'r')
    ax[4].plot(B_vals[int(num_samples/2):], M4_vals[int(num_samples/2):], 'b')
    ax[4].set_ylim(-1.0, 1.0)
    #ax[4].text(.3, -.5, "$K = .05$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))



    ax[4].set_ylabel("$M_{4}$", fontsize=13)
    ax[4].set_xlabel("K", fontsize=13)

    plt.show()

def plot_multi_data(sim_num_list, color_list, J_vals = []):
    fig = plt.figure(figsize=(7,5))
    ax0 = fig.add_subplot(1, 1, 1)

    for sim_num, col, j in zip(sim_num_list, color_list, J_vals):
        filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

        data = pickle.load(open(filename, 'rb'))
        num_samples = len(data[0])
        B_vals = data[0]
        M_vals = data[1]
        #M1_vals = data[2]
        #M2_vals = data[3]
        #M3_vals = data[4]

        ax0.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], color=col, linewidth=1.0, label=j)
        ax0.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], color =col, linewidth = 1.0)
        ax0.set_ylim(-1, 1)
        ax0.set_xlim(-.3,.3)
        #plt.setp(ax0.get_xticklabels(), visible=False)

        ax0.set_ylabel("$M$",fontsize=13)
    plt.legend(title = "$J_{AFM,1}$, $J_{AFM,2}$")
    ax0.set_xlabel("$B$", fontsize=13)


    plt.show()


def plot_data(sim_num):
    filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    #filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    data = pickle.load(open(filename, 'rb'))
    print(data)
    num_samples = len(data[0])
    B_vals = data[0]
    print("B vals: ", B_vals, "\n\n")
    M_vals = data[1]
    #M1_vals = data[2]
    #M2_vals = data[3]
    #M3_vals = data[4]
    #M4_vals = data[5]

    #fig, ax = plt.subplots(5, 1, figsize=(7,7), sharex=True)   # Magnetization per spin, all lattices
    #fig = plt.figure(figsize=(7,7))
    fig = plt.figure()
    #ax0 = fig.add_subplot(6, 1, (1,2))
    #ax0 = fig.add_subplot(4, 1, 1)
    ax0 =fig.add_subplot(1,1,1)
    #ax1 = fig.add_subplot(4, 1, 2, sharex = ax0)
    #ax2 = fig.add_subplot(4, 1, 3, sharex = ax0)
    #ax3 = fig.add_subplot(4, 1, 4, sharex = ax0)
    #ax4 = fig.add_subplot(6, 1, 6, sharex = ax0)



    #ax0.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    #ax0.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    ax0.plot(B_vals, M_vals)
    #ax0.set_ylim(-1, 1)
    ax0.set_xlim(0,5)
    #plt.setp(ax0.get_xticklabels(), visible=False)

    ax0.set_ylabel("$M$",fontsize=13)

    '''
    #plt.subplot(512)   # Magnetization per spin, first lattice
    ax1.plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ax1.plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    ax1.set_ylim(-1.0, 1.0)

    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1.set_ylabel("$M_{1}$",fontsize=13)
    #ax1.text(.11, -.5, "$K = 0.05$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    #ax1.text(.3, -.5, "$K = 0.05$\n$J_{AFM,1}=0.25$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    ax1.text(.22, -.5, "$K = 0.08$\n$J = 0.1$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))

    #plt.xticks([])
    #plt.xlabel("B")

    #plt.subplot(513)   # Magnetization per spin, second lattice
    ax2.plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ##ax[2].plot(B_vals[0:int(num_samples/4)], M2_vals[0:int(num_samples/4)], 'b')
    ##ax[2].plot(B_vals[int(num_samples/4):2*int(num_samples/4)], M2_vals[int(num_samples/4):2*int(num_samples/4)], 'r')

    ax2.plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    ##ax[2].plot(B_vals[2*int(num_samples/4):3*int(num_samples/4)], M2_vals[2*int(num_samples/4):3*int(num_samples/4)], 'b--')
    ##ax[2].plot(B_vals[3*int(num_samples/4):], M2_vals[3*int(num_samples/4):], 'r--')
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax2.set_ylim(-1.0, 1.0)


    ax2.set_ylabel("$M_{2}$",fontsize=13)
    #ax2.text(.11, -.5, "$K = 0.05$\n$J_{AFM,2} = 0.2$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    #ax2.text(.3, -.6, "$K = 0.05$\n$J_{AFM,2}=0.2$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    ax2.text(.22, -.5, "$K = 0.08$\n$J = 0.24$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))

    #plt.xticks([])

    #plt.xlabel("B")

    #plt.subplot(514)   # Magnetization per spin, third lattice

    ax3.plot(B_vals[0:int(num_samples/2)], M3_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ax3.plot(B_vals[int(num_samples/2):], M3_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    ##ax[3].plot(B_vals[0:int(num_samples/4)], M3_vals[0:int(num_samples/4)], 'b')
    ##ax[3].plot(B_vals[int(num_samples/4):2*int(num_samples/4)], M3_vals[int(num_samples/4):2*int(num_samples/4)], 'r')

        #ax[2].plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    ##ax[3].plot(B_vals[2*int(num_samples/4):3*int(num_samples/4)], M3_vals[2*int(num_samples/4):3*int(num_samples/4)], 'b--')
    ##ax[3].plot(B_vals[3*int(num_samples/4):], M3_vals[3*int(num_samples/4):], 'r--')
    plt.setp(ax3.get_xticklabels(), visible=False)

    ax3.set_ylim(-1.0, 1.0)


    ax3.set_ylabel("$M_{3}$",fontsize=13)
    #ax3.text(.11, 0, "$K = 0.05$\n$J_{AFM,3} = 0.2$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    #ax3.text(.3, -.6, "$K = 0.05$\n$J_{AFM,3}=0.2$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    ax3.text(.15, -.5, "$K = 0.08$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))

    #plt.xlabel("B")
    #plt.xticks([])
    '''

    '''
    #plt.subplot(515)   # Magnetization per spin, fourth lattice
    ax4.plot(B_vals[0:int(num_samples/2)], M4_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ax4.plot(B_vals[int(num_samples/2):], M4_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    ax4.set_ylim(-1.0, 1.0)


    ax4.text(.11, -.5, "$K = 0.05$", fontsize = 10, bbox=dict(edgecolor='black', fill = False, alpha=0.5))



    ax4.set_ylabel("$M_{4}$", fontsize=13)
    ax4.set_xlabel("$J_{AFM,1}$", fontsize=13)
    '''

    plt.show()

if __name__ == "__main__":
    #plot_paper_data()
    #plot_auto_data([179,180,181])
    #plot_data(155)
    visualize_lattice(397)

    #visualize_spins(231)
    #plot_multi_data([42,43,44,45,46,47,48,49], ['b','g','r','c','m','y','k', orange], [".15,.15",".18,.15",".21,.15",".24,.15",".1,.15",".1,.18",".1,.21",".1,.24"])
