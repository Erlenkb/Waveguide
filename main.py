import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D


###### Global Parameters ######
length = 2            #m
width = 0.2           #m
z_lim = 1
frequency = 2000      #Hz
c = 343               #m/s
rho = 1.25            #kgN/s

R_top    = 1  # Reflection factor top
R_bottom = 1  #
R_left   = 1  #
R_right  = 1  #

Resolution = 25   #Amount of steps per wavelength
################################




def _createWaveguideGrid():
    """
    Creates the number of gridpoints for the TLM matrix
    :return: The matrix as an N x M matrix
    """
    row = int(Resolution * width) - 1
    column = int(Resolution * length) - 1
    return np.zeros((row, column, 4)), np.zeros((row, column, 4))

def _create_scatter_node(p_matrix):
    """
    Creates the scatter matrix for a node
    :param p_matrix: Pressure array containing all incoming pressures for point i x j
    :return: Scattering array [4 elements]
    """
    S = np.zeros(4)
    s_matrix = np.multiply(p_matrix, np.ones((len(p_matrix), len(p_matrix))))
    for i in range(len(S)) : S[i] = np.sum(s_matrix[i])
    return S

def _create_harmonic(signal):
    return True

def _S_Matrix(pressure, scatter):
    """
    Creates the scatter matrix for the next timestep given the pressure for each point in the
    TLM matrix
    :param pressure: Pressure matrix
    :param scatter: Scatter matrix
    :return: Empty pressure matrix and the new scatter matrix
    """
    for i in range(pressure.shape[0]):
        for j in range(pressure.shape[1]): scatter[i][j] = _create_scatter_node(pressure[i][j])
    return pressure * 0, scatter

def _P_Matrix(scatter, pressure):
    """
    Creates the pressure for each node in the TLM matrix given the scatter matrix
    :param scatter: Scatter matrix
    :param pressure: Pressure matrix
    :return: The pressure matrix and the empty scatter matrix
    """
    list = [2,3,0,1]
    for i in range(scatter.shape[0]):
        for j in range(scatter.shape[1]):
            for p, s in enumerate(list):
                if (p == 0):
                    if (i == 0):
                        pressure[i][j][p] = R_top * scatter[i][j][p]
                    else: pressure[i][j][p] = scatter[i-1][j][s]

                elif (p == 1):
                    if (j == 0):
                        pressure[i][j][p] = R_left * scatter[i][j][p]
                    else: pressure[i][j][p] = scatter[i][j-1][s]

                elif (p == 2):
                    if (i == (scatter.shape[0] - 1)):
                        pressure[i][j][p] = R_bottom * scatter[i][j][p]
                    else: pressure[i][j][p] = scatter[i + 1][j][s]

                elif (p == 3):
                    if (j == (scatter.shape[1] - 1)):
                        pressure[i][j][p] = R_right * scatter[i][j][p]
                    else: pressure[i][j][p] = scatter[i][j + 1][s]
    return pressure, scatter*0

def _sum_pressurenode(pressure):
    """

    :param pressure: pressure array
    :return: Return the pressure for a given node
    """
    return .5 * np.sum(pressure)

def _heatmap(pressure):
    """
    Create a 2D heatmap for the TLM matrix to visualize the pressure value for each node
    :param pressure: Pressure matrix
    :return: 2D heatmap for the TLM matrix
    """
    heatmap = np.zeros((pressure.shape[0], pressure.shape[1]))
    for i in range(pressure.shape[0]):
        for j in range(pressure.shape[1]):
            heatmap[i][j] = _sum_pressurenode(pressure[i][j])
    return heatmap

def _harmonic_signal(time):
    return np.sin(2 * np.pi * frequency * time)

def _CreateAnimation(source, num_it, filename):
    P_grid, S_grid = _createWaveguideGrid()

    timestep = round((1 / (Resolution * frequency)), 7)
    if (source == "Point") : P_grid[0][0] = np.array([.25,.25,.25,.25])
    filename_early = []
    filename_late = []

    for i in range(num_it):
        sig = round(0.5 * _harmonic_signal(i * timestep),4)
        if (source == "Harmonic") :
            P_grid[0][0][0] = 0.5 * sig #+ P_grid[0][0][0]
            P_grid[0][0][3] = 0.5 * sig #+ P_grid[0][0][3]
        P_grid, S_grid = _S_Matrix(P_grid, S_grid)

        P_grid, S_grid = _P_Matrix(S_grid, P_grid)
        print ("It. number: {0}".format(i))

        if (i > 100 and i < 200):
            Z = _heatmap(P_grid)
            filename_early.append("{0}.png".format(i))

        elif (i > 2000 and i < 2200) :
            Z = _heatmap(P_grid)
            filename_late.append("{0}.png".format(i))

        else : continue

        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection="3d")

        x = np.linspace(0, length, P_grid.shape[1])
        y = np.linspace(0, width, P_grid.shape[0])
        X, Y = np.meshgrid(x,y)

        ax.set_zlim(-z_lim, z_lim)
        ax.set_ylim(0, width)
        ax.set_xlim(0, length)

        ax.set_xlabel("length [m]")
        ax.set_ylabel("width [m]")
        ax.set_zlabel("Amplitude")

        x_scale = 2
        y_scale = 0.4
        z_scale = 1

        scale = np.diag([x_scale, y_scale, z_scale, 1])
        scale = scale * (1 / scale.max())

        scale[3, 3] = 1

        def short_proj():
            return np.dot(Axes3D.get_proj(ax), scale)

        ax.get_proj = short_proj

        ax.plot_surface(X, Y, _heatmap(P_grid))
        plt.savefig("{0}".format(i))

        plt.close(fig)

    return filename_early, filename_late







if __name__ == '__main__':
    _CreateAnimation("Point", 400, "cool")



