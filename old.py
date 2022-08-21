import pandas as pd
import os
import seaborn as sb
import numpy as np
import imageio
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import scipy.signal as signal

######## IMPORTANT PARAMETERS ########

TLM_Row_size = 0.2  # In meters
TLM_Column_size = 2  # In meters
frequency = 500  # Hz
# Timesteps =    #    TLM_Column_size * 100
R_top = 1        #
R_bottom = 1     #
R_left = 1       #
R_right = 1      #
Resolution = 25  # Defines the amount of steps per wavelength


def TLM_SIZE():
    """
    Creates the number of gridpoints for the TLM matrix
    :return: The matrix as row x column
    """
    wave = 343 / frequency
    stepsize = wave / Resolution
    step_row = TLM_Row_size / stepsize
    step_column = TLM_Column_size / stepsize
    return int(step_row), int(step_column)


def TLM(row, column):
    """
    Creates a 3D matrix with the size of the TLM Matrix
    :param row: Row size of the TLM matrix
    :param column: Column size of the TLM matrix
    :return: returns 2 3D matrices with the TLM diameters
    """
    return np.zeros((row, column, 4)), np.zeros((row, column, 4))


def create_harmonic(frequency, k):
    """
    Creates the harmonic signal
    :param frequency: Source signal frequency
    :param k: k is the timestep count for the signal
    :return: returns the signal value as a single scalar value.
    """
    wave = 343 / frequency
    stepsize = wave / Resolution
    step_length = k * (stepsize / 343)

    signal = np.sin(2 * np.pi * frequency * step_length)
    signal_array = np.zeros((Resolution, 4))
    # for i in range(signal_array.shape[0]):
    #    signal_array[i] = np.array([signal[i],0,0,0])
    return signal


def scatter_node(pressure):
    """
    Creates the scatter matrix for a node
    :param pressure: Pressure array containing all incoming pressures for a point i,j
    :return: scattering array [4 elements]
    """
    S = np.zeros(4)
    s_matrix = np.multiply(pressure, np.ones((len(pressure), len(pressure))) + np.eye(len(pressure)) * (-2)) * 0.5
    for i in range(len(S)):
        S[i] = np.sum(s_matrix[i])
    return S


def create_S_grid(Pressure, Scatter):
    """
    Creates the scatter matrix for the next timestep given the pressure for each point in the TLM matrix
    :param Pressure: Pressure matrix
    :param Scatter: Scatter matrix
    :return: empty pressure matrix and the scatter matrix
    """
    for i in range(Pressure.shape[0]):
        for j in range(Pressure.shape[1]):
            Scatter[i][j] = scatter_node(Pressure[i][j])
    return Pressure * 0, Scatter


def create_p_grid(Scatter, Pressure):
    """
    Creates the pressure for each node in the TLM matrix given the scatter
    :param Scatter: Scatter matrix
    :param Pressure: Pressure matrix
    :return: The pressure and empty scatter matrix
    """

    list = [2, 3, 0, 1]
    for i in range(Scatter.shape[0]):
        for j in range(Scatter.shape[1]):
            for p, s in enumerate(list):
                if (p == 0):
                    if (i == 0):
                        Pressure[i][j][p] = R_top * Scatter[i][j][p]
                    else:
                        Pressure[i][j][p] = Scatter[i - 1][j][s]
                elif (p == 1):
                    if (j == 0):
                        Pressure[i][j][p] = R_left * Scatter[i][j][p]
                    else:
                        Pressure[i][j][p] = Scatter[i][j - 1][s]
                elif (p == 2):
                    if (i == (Scatter.shape[0] - 1)):
                        Pressure[i][j][p] = R_bottom * Scatter[i][j][p]
                    else:
                        Pressure[i][j][p] = Scatter[i + 1][j][s]
                elif (p == 3):
                    if (j == (Scatter.shape[1] - 1)):
                        Pressure[i][j][p] = R_right * Scatter[i][j][p]
                    else:
                        Pressure[i][j][p] = Scatter[i][j + 1][s]

    return Pressure, Scatter * 0


def sum_pressure(pressure):
    """

    :param pressure: Pressure array
    :return: return the total pressure at spesific node
    """
    return 0.5 * np.sum(pressure)


def heatmap(Pressure):
    """
    Creates a 2D heatmap for the TLM matrix to visualize the pressure value for each node.
    :param Pressure: Pressure matrix
    :return: Returns a 2D heatmap for the TLM matrix
    """
    Heat = np.zeros((Pressure.shape[0], Pressure.shape[1]))
    for i in range(Pressure.shape[0]):
        for j in range(Pressure.shape[1]):
            Heat[i][j] = sum_pressure(Pressure[i][j])
    return Heat


def iterate_matrix(num, c=True):
    """
    Creates the TLM matrix and iterate the matrix for a given amount of timesteps.
    It also creates num amount of pictures for the matrix and creates a Gif file for the simulation. Some changes where made to the code due to time restrictions, but the main
    Objective for this code is to run all the functions, create the matrix and iterate the source signal and watch its wave propogate through the waveguide.

    It did print the TLM matrix both in 3D and 2D, but focused at the end for 3D, hance the bit messy code.
    :param num: Number of iterations for the simulation
    :param c: if c==True it will create a outgoing wave in all branches for the lower-left corner node.
    :return:
    """
    row, column = TLM_SIZE()
    P_grid, S_grid = TLM(row, column)

    if c:
        P_grid[0][0] = np.array([0.25, 0.25, 0.25, 0.25])

    counter = 0
    filenames = []
    fileCount = 0
    fft_array = []
    for i in range(num):
        sig = create_harmonic(frequency, i)
        if (c == False):
            P_grid[0][0][0] = 0.5 * create_harmonic(frequency, i) + P_grid[0][0][0]
            P_grid[0][0][3] = 0.5 * create_harmonic(frequency, i) + P_grid[0][0][3]
        P_grid, S_grid = create_S_grid(P_grid, S_grid)

        P_grid, S_grid = create_p_grid(S_grid, P_grid)

        counter += 1



        """ 
        if(counter>=9): counter=0
        if (fileCount==3):

            fig = plt.figure(figsize=(10,1))
            sb.heatmap(heatmap(P_grid),-0.6,0.6)

            plt.savefig("{0}.png".format(i))
            filenames.append("{0}.png".format(i))
            plt.close(fig)
            fileCount=0
        """

        if (i > 0):
            Z = heatmap(P_grid)
            fft_array.append(Z[1][P_grid.shape[1] - 1])

        if (i > 0):
            fig = plt.figure(figsize=(10, 10))
            ax = fig.gca(projection='3d')

            Z = heatmap(P_grid)
            x = np.linspace(0, 2, P_grid.shape[1])
            y = np.linspace(0, 0.2, P_grid.shape[0])
            X, Y = np.meshgrid(x, y)
            ax.set_zlim(-0.3, 0.3)
            ax.set_ylim(0, 0.2)
            ax.set_xlim(0, 2)

            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            ax.set_zlabel("Amplitude [P]")

            x_scale = 2
            y_scale = 0.4
            z_scale = 1

            scale = np.diag([x_scale, y_scale, z_scale, 1.0])
            scale = scale * (1.0 / scale.max())
            scale[3, 3] = 1.0

            def short_proj():
                return np.dot(Axes3D.get_proj(ax), scale)

            ax.get_proj = short_proj

            ax.plot_surface(X, Y, heatmap(P_grid))
            plt.savefig("{0}.png".format(i))
            filenames.append("{0}.png".format(i))
            # plt.show()
            plt.close(fig)
            fileCount += 1
    wave = 343 / frequency
    stepsize = wave / Resolution
    step_length = (stepsize / 343)
    sp = np.fft.fft(fft_array, 44100)
    sp = np.trim_zeros(sp, trim="fb")
    freq = np.fft.fftfreq(n=len(sp), d=step_length)
    plt.plot(freq, np.abs(sp), color="dimgray", label="FFT of recorded signal")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.grid()
    plt.xscale("log")
    plt.legend(fontsize="x-large")
    peaks = signal.find_peaks(sp)
    print(peaks)
    plt.xlim(150, 2000)

    plt.show()

    with imageio.get_writer("{0}Hz{1}Res R{2}.gif".format(frequency, Resolution, R_right), mode="I") as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
    print("Gif created\n")
    for filename in set(filenames):
        os.remove(filename)
    print("done")


if __name__ == '__main__':
    # P_init = create_harmonic(2000,10)
    Pressure_array = [1, 0, 0, 0]

    iterate_matrix(150, c=False)



