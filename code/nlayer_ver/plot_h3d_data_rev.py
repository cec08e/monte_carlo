import pickle
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import Normalize

from mpl_toolkits.mplot3d import Axes3D
from itertools import cycle
from functools import reduce
import numpy as np
import csv
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
[0.35, -0.417266],
[0.4,-0.336577],
[0.45,-0.26326],
[0.5,-0.188524],
[0.55,-0.137358],
[0.6,-0.089943],
[0.65,-0.062394],
[0.7,-0.065756],
[0.75,-0.058107],
[0.8,-0.057072],
[0.85,-0.070517],
[0.9,-0.079535],
[0.95,-0.079105],
[1,-0.07509],
[1.05,-0.078691],
[1.1,-0.07878],
[1.15,-0.070854],
[1.2,-0.072561],
[1.25,-0.074923],
[1.3,-0.064094],
[1.35,-0.064082],
[1.4,-0.060505],
[1.45,-0.0604],
[1.5,-0.05758],
[1.55,-0.050159],
[1.6,-0.054178],
[1.65,-0.052888],
[1.7,-0.042017],
[1.75,-0.040414],
[1.8,-0.037921],
[1.85,-0.037285],
[1.9,-0.036258],
[1.95,-0.034623],
[2,-0.030835],
[2.05,-0.034204],
[2.1,-0.033895],
[2.15,-0.025037],
[2.2,-0.027537],
[2.25,-0.02288],
[2.3,-0.028378],
[2.35,-0.027442],
[2.4,-0.021257],
[2.45,-0.019514]
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
[0.25,-0.192012],
[0.3,-0.216405],
[0.35,-0.229675],
[0.4,-0.222256],
[0.45,-0.200866],
[0.5,-0.14635],
[0.55,-0.106145],
[0.6,-0.073286],
[0.65,-0.050359],
[0.7,-0.056516],
[0.75,-0.048562],
[0.8,-0.063556],
[0.85,-0.069976],
[0.9,-0.071689],
[0.95,-0.083587],
[1,-0.085832],
[1.05,-0.084987],
[1.1,-0.089011],
[1.15,-0.086372],
[1.2,-0.085122],
[1.25,-0.075349],
[1.3,-0.077443],
[1.35,-0.068196],
[1.4,-0.063161],
[1.45,-0.059969],
[1.5,-0.061689],
[1.55,-0.054922],
[1.6,-0.057585],
[1.65,-0.053678],
[1.7,-0.048427],
[1.75,-0.051302],
[1.8,-0.049519],
[1.85,-0.046929],
[1.9,-0.039484],
[1.95,-0.034636],
[2,-0.03747],
[2.05,-0.034606],
[2.1,-0.026392],
[2.15,-0.03408],
[2.2,-0.03187],
[2.25,-0.025317],
[2.3,-0.026577],
[2.35,-0.023657],
[2.4,-0.019719],
[2.45,-0.021208],
[2.5,-0.029074]
]
'''
[2.55,-0.017023],
[2.6,-0.023096],
[2.65,-0.019355],
[2.7,-0.011201],
[2.75,-0.0185],
[2.8,-0.017305],
[2.85,-0.018676],
[2.9,-0.004097],
[2.95,-0.01412],
[3,-0.011721],
[3.05,-0.017301],
[3.1,-0.008103],
[3.15,-0.021561],
[3.2,-0.008571],
[3.25,-0.01452],
[3.3,-0.009196],
[3.35,-0.007861],
[3.4,-0.011173],
[3.45,-0.006933],
[3.5,-0.000385],
[3.55,-0.008365],
[3.6,-0.011643],
[3.65,-0.003953],
[3.7,-0.011993],
[3.75,-0.011851],
[3.8,-0.012194],
[3.85,-0.00402],
[3.9,-0.005391],
[3.95,-0.009874],
[4,-0.008639]
]
'''
# 1.15 #
Q_v_B_AF_6 = [
[0.025,-0.827989],
[0.05,-1.073296],
[0.075,-1.102854],
[0.1,-1.071698],
[0.125,-1.01304],
[0.15,-0.955409],
[0.175,-0.866776],
[0.2,-0.802581],
[0.225,-0.74446],
[0.25,-0.680952],
[0.5,-0.21809],
[0.7,-0.060397],
[1,-0.070624],
[1.5,-0.05254],
[2,-0.0289],
[2.5,-0.020938]
]
#1.0#
Q_v_B_AF_7 = [
[0.025,-1.09776],
[0.05,-1.264299],
[0.075,-1.242208],
[0.1,-1.161864],
[0.125,-1.079426],
[0.15,-0.962844],
[0.175,-0.892005],
[0.2,-0.813198],
[0.225,-0.743337],
[0.25,-0.686463],
[0.5,-0.224855],
[0.7,-0.062522],
[1,-0.064029],
[1.5,-0.052785],
[2,-0.022777],
[2.5,-0.01776]
]
#.85#
Q_v_B_AF_8 = [
[0.025,-1.249953],
[0.05,-1.319897],
[0.075,-1.249659],
[0.1,-1.144456],
[0.125,-1.026558],
[0.15,-0.93536],
[0.175,-0.844363],
[0.2,-0.767609],
[0.225,-0.6874],
[0.25,-0.634629],
[0.5,-0.231686],
[0.7,-0.06914],
[1,-0.058064],
[1.5,-0.032891],
[2,-0.022592],
[2.5,-0.0041]
]

#.6#
Q_v_B_AF_9 = [
[0.025,-1.171289],
[0.05,-1.161554],
[0.075,-1.045476],
[0.1,-0.930633],
[0.125,-0.825562],
[0.15,-0.734404],
[0.175,-0.661428],
[0.2,-0.596549],
[0.225,-0.53888],
[0.25,-0.497868],
[0.5,-0.180442],
[0.7,-0.057785],
[1,-0.042142],
[1.5,-0.026648],
[2,-0.012311],
[2.5,-0.006715]
]

#.45#
Q_v_B_AF_10 = [
[0.025,-0.956847],
[0.05,-0.936139],
[0.075,-0.83794],
[0.1,-0.726709],
[0.125,-0.647931],
[0.15,-0.575262],
[0.175,-0.512059],
[0.2,-0.461659],
[0.225,-0.413983],
[0.25,-0.37934],
[0.5,-0.142451],
[0.7,-0.040582],
[1,-0.033778],
[1.5,-0.022514],
[2,-0.012181],
[2.5,-0.016843]
]

#.3#
Q_v_B_AF_11 = [
[0.05,-0.654178],
[0.5,-0.09888],
[0.7,-0.03247],
[1,-0.023327],
[1.5,-0.013573],
[2,-0.010035],
[2.5,-0.004772]
]

#.15#
Q_v_B_AF_12 = [
[0.05,-0.33393],
[0.5,-0.055594],
[0.7,-0.017979],
[1,-0.000882],
[1.5,0.004113],
[2,-0.006945],
[2.5,-0.001231]
]

#0#
Q_v_B_AF_13 = [
[0.05,-0.004269],
[0.5,0.005828],
[0.7,0.006858],
[1,0.007421],
[1.5,0.007936],
[2,0.004767],
[2.5,-0.001819]
]

def plot_heat_map():
    filename = "AF_results/AF_data.csv"
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        data = []
        #print(reader.fieldnames)
        #for fn in reader.fieldnames:
        #    data.append([])


        for i,row in enumerate(reader):
            if i == 0:
                continue
            if float(row[0]) > 1.0:
                break
            #for i,n in enumerate(reader.fieldnames):
            #    data[i].append(row[n])

            data.append(row[1:])

        print("data: ", data[:][2])

        fig, ax = plt.subplots(1, 1, figsize=(10,10))   # Magnetization per spin, all lattices

        ax.imshow([[float(item) for item in np.array(data)[:,i]] for i in range(len(data[0]))], cmap=plt.get_cmap('hot_r'), interpolation='bilinear')
        ax.set_xlabel("T")
        ax.set_ylabel("B")
        plt.show()



def visualize_lattice(sim_num):
    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

    data = pickle.load(open(filename, 'rb'))
    print(data)
    plt.imshow([[item[2] for item in row] for row in data[0]], cmap=plt.get_cmap('Spectral'))
    plt.show()

def visualize_bilayer_lattice(sim_num):
    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    filename = "periodic_inter/sim_results/sim_" + str(sim_num) + ".pickle"

    data = pickle.load(open(filename, 'rb'))
    f, ax = plt.subplots(2)
    norm = Normalize(vmin=-1.0, vmax = 1.0)
    print(data)
    for i, layer in enumerate(data):
        ax[i].imshow([[item[2] for item in row] for row in layer], norm = norm, cmap=plt.get_cmap('Spectral'))
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

def visualize_TC_density(sim_num):
    filename = "dirac_results/sim_" + str(sim_num) + ".pickle"

    data = pickle.load(open(filename, 'rb'))
    print(data)
    # triangulate lattice and calculate the TC vals and plot
    TC_points = []
    x = []
    y = []
    t = []
    for i, row in enumerate(data[0]):
        for j,spin in enumerate(row):
            print("i,j: ", i, ", ", j)
            #index of point given by i*(len(row)) + j
            x.append(i)
            y.append(j)
            # Calc TC for (i,j) --> (i, j-1) ---> (i+1, j) and plot at (i+1/3, j-1/3)
            if j-1 >= 0 and i+1 < len(data[0]):
                #TC_points.append((i+(1/3.0), j-(1/3.0), calc_SA(data[0][i][j], data[0][i][(j-1)%len(row)], data[0][(i+1)%len(data[0])][j])))
                TC_points.append(calc_SA(data[0][i][j], data[0][i][(j-1)%len(row)], data[0][(i+1)%len(data[0])][j]))
                t.append((i*len(row) + j, i*len(row) + j - 1, (i+1)*len(row) + j))
                print("Appending triangle: ", (i*len(row) + j, i*len(row) + j - 1, (i+1)*len(row) + j))

            # Calc TC for (i,j) --> (i, j+1) ---> (i-1, j) and plot at (i-1/3, j+1/3)
            if i-1 >= 0 and j+1 < len(row):
                #TC_points.append((i-(1/3.0), j+(1/3.0), calc_SA(data[0][i][j], data[0][i][(j+1)%len(row)], data[0][(i-1)%len(data[0])][j])))
                TC_points.append(calc_SA(data[0][i][j], data[0][i][(j+1)%len(row)], data[0][(i-1)%len(data[0])][j]))
                t.append((i*len(row) + j, i*len(row) + j + 1, (i-1)*len(row) + j))

    print(len(t))
    print(len(TC_points))
    print(reduce(lambda x,y: x+y, TC_points))

    triang = Triangulation(x,y,t)
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    tpc = ax.tripcolor(triang, TC_points, cmap=plt.get_cmap('PiYG'))
    fig.colorbar(tpc)
    plt.gca().invert_yaxis()

    plt.show()

    '''
    cmap = matplotlib.cm.get_cmap('cool')
    for elem in TC_points:
        plt.plot(elem[0], elem[1], cmap=plt.get_cmap('Spectral'))
    plt.show()
    '''

def visualize_bilayer_TC_density(sim_num):
    filename = "periodic_inter/sim_results/sim_" + str(sim_num) + ".pickle"

    data = pickle.load(open(filename, 'rb'))
    print(data)
    # triangulate lattice and calculate the TC vals and plot
    TC_points = []
    x = []
    y = []
    t = []
    fig, ax = plt.subplots(2)
    TC_total = 0.0

    for k, layer in enumerate(data):
        for i, row in enumerate(layer):
            for j,spin in enumerate(row):
                #print("i,j: ", i, ", ", j)
                #index of point given by i*(len(row)) + j
                x.append(i)
                y.append(j)

                # Calc TC for (i,j) --> (i, j-1) ---> (i+1, j) and plot at (i+1/3, j-1/3)
                TC_total += calc_SA(data[k][i][j], data[k][i][(j-1)%len(row)], data[k][(i+1)%len(data[0])][j])
                if j-1 >= 0 and i+1 < len(data[k]):
                    #TC_points.append((i+(1/3.0), j-(1/3.0), calc_SA(data[0][i][j], data[0][i][(j-1)%len(row)], data[0][(i+1)%len(data[0])][j])))
                    TC_points.append(calc_SA(data[k][i][j], data[k][i][(j-1)%len(row)], data[k][(i+1)%len(data[0])][j]))
                    t.append((i*len(row) + j, i*len(row) + j - 1, (i+1)*len(row) + j))
                    #print("Appending triangle: ", (i*len(row) + j, i*len(row) + j - 1, (i+1)*len(row) + j))

                # Calc TC for (i,j) --> (i, j+1) ---> (i-1, j) and plot at (i-1/3, j+1/3)
                TC_total += calc_SA(data[k][i][j], data[k][i][(j+1)%len(row)], data[k][(i-1)%len(data[0])][j])
                if i-1 >= 0 and j+1 < len(row):
                    #TC_points.append((i-(1/3.0), j+(1/3.0), calc_SA(data[0][i][j], data[0][i][(j+1)%len(row)], data[0][(i-1)%len(data[0])][j])))
                    TC_points.append(calc_SA(data[k][i][j], data[k][i][(j+1)%len(row)], data[k][(i-1)%len(data[0])][j]))
                    t.append((i*len(row) + j, i*len(row) + j + 1, (i-1)*len(row) + j))

        #print(len(t))
        #print(len(TC_points))
        print("TC: ", TC_total)
        TC_total = 0.0

        triang = Triangulation(x,y,t)
        ax[k].set_aspect('equal')
        print(TC_points)
        norm = Normalize( vmin = -1*(max([-1*min(TC_points), max(TC_points)])), vmax = max([-1*min(TC_points), max(TC_points)]))

        tpc = ax[k].tripcolor(triang, TC_points, cmap=plt.get_cmap('PiYG'), norm=norm)
        fig.colorbar(tpc, ax = ax[k])
        plt.gca().invert_yaxis()
        TC_points = []
        x = []
        y = []
        t = []

    #fig.colorbar(tpc)


    plt.show()

def visualize_bilayer_quiver(sim_num):
    fig, ax = plt.subplots()
    filename = "periodic_inter/sim_results/sim_" + str(sim_num) + ".pickle"
    data = pickle.load(open(filename, 'rb'))
    U = [[item[0] for item in row] for row in data[0]]
    V = [[item[1] for item in row] for row in data[0]]
    q = ax.quiver( U, V)
    plt.gca().invert_yaxis()

    plt.show()

def calc_SA(n1, n2, n3):
    n2_cross_n3 = np.cross(n2, n3)
    n1_dot_n2 = np.dot(n1, n2)
    n2_dot_n3 = np.dot(n2, n3)
    n3_dot_n1 = np.dot(n3, n1)


    rho = np.sqrt(2*(1+n1_dot_n2)*(1+n2_dot_n3)*(1+n3_dot_n1));
    temp = complex((1.0/rho)*(1 + n1_dot_n2 + n2_dot_n3 + n3_dot_n1), (1/rho)*np.dot(n1, n2_cross_n3))

    Omega = (2*np.log(temp).imag)/(4*np.pi)
    #print("Omega: ", Omega)

    return Omega



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
    AF_line_4, = plt.plot([x[0] for x in Q_v_B_AF_4], [x[1]*(1000/400.0) for x in Q_v_B_AF_4], 'og', label="$D=.3$, $B=1.3$")
    AF_line_5, = plt.plot([x[0] for x in Q_v_B_AF_5], [x[1]*(1000/400.0) for x in Q_v_B_AF_5], 'ob', label="$D=.3$, $B=1.45$")
    AF_line_6, = plt.plot([x[0] for x in Q_v_B_AF_6], [x[1]*(1000/400.0) for x in Q_v_B_AF_6], 'or', label="$D=.3$, $B=1.15$")
    AF_line_7, = plt.plot([x[0] for x in Q_v_B_AF_7], [x[1]*(1000/400.0) for x in Q_v_B_AF_7], 'om', label="$D=.3$, $B=1.0$")
    AF_line_8, = plt.plot([x[0] for x in Q_v_B_AF_8], [x[1]*(1000/400.0) for x in Q_v_B_AF_8], 'oc', label="$D=.3$, $B=.85$")
    AF_line_9, = plt.plot([x[0] for x in Q_v_B_AF_9], [x[1]*(1000/400.0) for x in Q_v_B_AF_9], 'ok', label="$D=.3$, $B=.6$")
    AF_line_10, = plt.plot([x[0] for x in Q_v_B_AF_10], [x[1]*(1000/400.0) for x in Q_v_B_AF_10], 'oy', label="$D=.3$, $B=.45$")
    AF_line_11, = plt.plot([x[0] for x in Q_v_B_AF_11], [x[1]*(1000/400.0) for x in Q_v_B_AF_11], '^g', label="$D=.3$, $B=.3$")
    AF_line_12, = plt.plot([x[0] for x in Q_v_B_AF_12], [x[1]*(1000/400.0) for x in Q_v_B_AF_12], '^b', label="$D=.3$, $B=.15$")
    AF_line_13, = plt.plot([x[0] for x in Q_v_B_AF_13], [x[1]*(1000/400.0) for x in Q_v_B_AF_13], '^r', label="$D=.3$, $B=.0$")
    plt.legend(handles=[AF_line_13, AF_line_12, AF_line_11,AF_line_10, AF_line_9, AF_line_8, AF_line_7, AF_line_6, AF_line_4, AF_line_5], title="N = 20")


    plt.show()

def plot_AF_data():
    filename = "AF_results/AF_data.csv"
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        data = []
        print(reader.fieldnames)
        for fn in reader.fieldnames:
            data.append([])

        for row in reader:
            if float(row[reader.fieldnames[0]]) > 1.0:
                break
            for i,n in enumerate(reader.fieldnames):
                data[i].append(row[n])

        print(data)

        plt.ylabel("Q/1000 spins")
        plt.xlabel("T")
        color_cycle = cycle('rgbkmcy')
        handle_list = []
        for i in range(1,len(reader.fieldnames)):
            temp, =  plt.plot([float(x) for x in data[0]], [float(x)*(1000/400.0) for x in data[i]], '^'+next(color_cycle), label=reader.fieldnames[i])
            handle_list.append(temp)
    plt.legend(handles=handle_list, title="N = 20")


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
    #filename = "dirac_results/sim_" + str(sim_num) + ".pickle"
    filename = "M_R_results/sim_" + str(sim_num) + ".pickle"


    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    data = pickle.load(open(filename, 'rb'))
    print(data)
    num_samples = len(data[0])
    B_vals = data[0]
    print("B vals: ", B_vals, "\n\n")
    M_vals = data[1]
    #M1_vals = data[2]
    TC_vals = [ x*(1000/400.0) for x in data[3]]
    #M2_vals = data[3]
    #M3_vals = data[4]
    #M4_vals = data[5]

    #fig, ax = plt.subplots(5, 1, figsize=(7,7), sharex=True)   # Magnetization per spin, all lattices
    #fig = plt.figure(figsize=(7,7))
    fig = plt.figure()
    #ax0 = fig.add_subplot(6, 1, (1,2))
    #ax0 = fig.add_subplot(4, 1, 1)
    ax0 = fig.add_subplot(2,1,1)
    ax1 = fig.add_subplot(2, 1, 2, sharex = ax0)
    #ax2 = fig.add_subplot(4, 1, 3, sharex = ax0)
    #ax3 = fig.add_subplot(4, 1, 4, sharex = ax0)
    #ax4 = fig.add_subplot(6, 1, 6, sharex = ax0)



    ax0.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ax0.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    #ax0.plot(B_vals, M_vals)
    ax0.set_ylim(-1, 1)
    ax0.set_xlim(-.2,.2)
    #plt.setp(ax0.get_xticklabels(), visible=False)

    ax0.set_ylabel("$M$",fontsize=13)

    ax1.plot(B_vals[0:int(num_samples/2)], TC_vals[0:int(num_samples/2)], color=orange, linewidth=2.0)
    ax1.plot(B_vals[int(num_samples/2):], TC_vals[int(num_samples/2):], color =olive, linewidth = 2.0)
    #ax1.set_ylim(-1.0, 1.0)

    #plt.setp(ax1.get_xticklabels(), visible=False)

    ax1.set_ylabel("$Q/1000 spins$",fontsize=13)

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
    #plot_data(10)
    #plot_heat_map()
    #plot_AF_data()
    visualize_bilayer_lattice(268)
    #visualize_lattice(469)
    #visualize_TC_density(469)
    visualize_bilayer_TC_density(268)
    visualize_bilayer_quiver(268)


    #visualize_spins(231)
    #plot_multi_data([42,43,44,45,46,47,48,49], ['b','g','r','c','m','y','k', orange], [".15,.15",".18,.15",".21,.15",".24,.15",".1,.15",".1,.18",".1,.21",".1,.24"])
