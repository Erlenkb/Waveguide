import numpy as np
import matplotlib.pyplot as plt

Resolution = 25
c = 343
frequency = 2000

timestep = round((1 / Resolution) / c, 7)


time = np.linspace(0,1, Resolution * frequency)

val = np.sin(2 * np.pi * frequency * time)

plt.plot(time, val)
plt.show()

print('Hello world')