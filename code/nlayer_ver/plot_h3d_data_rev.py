import pickle
from matplotlib import pyplot as plt
#import matplotlib.gridspec as gridspec

orange = "#FF8000"
olive = "#008000"

def visualize_lattice(sim_num):
    filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    data = pickle.load(open(filename, 'rb'))
    print(data)
    plt.imshow([[item[2] for item in row] for row in data[0]], cmap=plt.get_cmap('Spectral'))
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
    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    data = pickle.load(open(filename, 'rb'))
    print(data)
    num_samples = len(data[0])
    B_vals = data[0]
    print("B vals: ", B_vals, "\n\n")
    M_vals = data[1]
    M1_vals = data[2]
    M2_vals = data[3]
    M3_vals = data[4]
    #M4_vals = data[5]

    #fig, ax = plt.subplots(5, 1, figsize=(7,7), sharex=True)   # Magnetization per spin, all lattices
    fig = plt.figure(figsize=(7,7))
    #ax0 = fig.add_subplot(6, 1, (1,2))
    ax0 = fig.add_subplot(4, 1, 1)
    ax1 = fig.add_subplot(4, 1, 2, sharex = ax0)
    ax2 = fig.add_subplot(4, 1, 3, sharex = ax0)
    ax3 = fig.add_subplot(4, 1, 4, sharex = ax0)
    #ax4 = fig.add_subplot(6, 1, 6, sharex = ax0)



    ax0.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ax0.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    ax0.set_ylim(-1, 1)
    ax0.set_xlim(-.3,.3)
    #plt.setp(ax0.get_xticklabels(), visible=False)

    ax0.set_ylabel("$M$",fontsize=13)


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
    #plot_data(50)
    visualize_lattice(63)
    #plot_multi_data([42,43,44,45,46,47,48,49], ['b','g','r','c','m','y','k', orange], [".15,.15",".18,.15",".21,.15",".24,.15",".1,.15",".1,.18",".1,.21",".1,.24"])
