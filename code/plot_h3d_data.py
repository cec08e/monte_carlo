import pickle
from matplotlib import pyplot as plt


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

    fig, ax = plt.subplots(5, 1, figsize=(8,10), sharex=True)   # Magnetization per spin, all lattices
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



def plot_data(sim_num):
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

    fig, ax = plt.subplots(5, 1, figsize=(8,10), sharex=True)   # Magnetization per spin, all lattices
    ax[0].plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    ax[0].plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')
    ax[0].set_ylim(-1.0, 1.0)
    ax[0].set_xlim(-.375,.425)

    ax[0].set_ylabel("M",fontsize=13)
    #plt.xlabel("B")
    #plt.xticks([])

    #plt.subplot(512)   # Magnetization per spin, first lattice
    ax[1].plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    ax[1].plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')
    ax[1].set_ylim(-1.0, 1.0)


    ax[1].set_ylabel("$M_{1}$",fontsize=13)
    ax[1].text(.3, -.5, "$K = .2$\n$J_{AFM} = .3$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))
    #plt.xticks([])
    #plt.xlabel("B")

    #plt.subplot(513)   # Magnetization per spin, second lattice
    ax[2].plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    ax[2].plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    ax[2].set_ylim(-1.0, 1.0)


    ax[2].set_ylabel("$M_{2}$",fontsize=13)
    ax[2].text(.3, -.5, "$K = .05$\n$J_{AFM} = .05$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))

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
    ax[4].text(.3, -.5, "$K = .05$", fontsize = 11, bbox=dict(edgecolor='black', fill = False, alpha=0.5))



    ax[4].set_ylabel("$M_{4}$", fontsize=13)
    ax[4].set_xlabel("B", fontsize=13)

    plt.show()





if __name__ == "__main__":
    plot_k_data(152)
