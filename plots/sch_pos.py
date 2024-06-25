import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("sch_point.txt")

x = data[:,0]

rew1 = data[:,1]
imw1 = data[:,2]

plt.plot(x, rew1)
plt.plot(x, imw1)

plt.show()