import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def f_exp(x,a,b):

    return a*np.exp(b*x)

xAIM = np.linspace(1, 10, 10)
xIAIM = np.linspace(1, 21, 21)

x = np.linspace(1,22,100)

dataAIM = np.loadtxt("pt_time_aim.txt", skiprows=1)
dataIAIM = np.loadtxt("pt_time_iaim.txt", skiprows=1)

poptAIM, pcovAIM = curve_fit(f_exp, dataAIM[:,0], dataAIM[:,1])
poptIAIM, pcovIAIM = curve_fit(f_exp, dataIAIM[:,0], dataIAIM[:,1])

fig = plt.figure()
ax = fig.gca()

# INTEGER AXIS TICKS
ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
# LOGARITHMIC SCALE ON Y AXIS
#ax.set_yscale("log")
# AXIS TITLES
ax.set_xlabel(r"Iteraciones $N$")
ax.set_ylabel(r"Tiempo de cálculo $t$ [s]")
# LOGARITHMIC SCALE ON Y AXIS
ax.set_yscale("log")
# AXIS LIMITS
ax.set_xlim(0.5, 22.5)
ax.set_ylim(10e-2,50e0)
# PLOT TITLE
ax.set_title("Complejidad temporal de los métodos AIM/IAIM")

ax.set_xticks([2,4,6,8,10,12,14,16,18,20,22])

ax.plot(x, f_exp(x, *poptAIM), alpha=0.55, color="tab:blue", linestyle="-")
ax.plot(x, f_exp(x, *poptIAIM), alpha=0.55, color="tab:orange", linestyle="-")

ax.plot(dataAIM[:,0], dataAIM[:,1], marker=".", linestyle="", label="Método AIM", color="tab:blue")
ax.plot(dataIAIM[:,0], dataIAIM[:,1], marker=".", linestyle="", label="Método IAIM", color="tab:orange")

#ax.plot(xAIM, f_exp(xAIM, *poptAIM), marker="", linestyle="-")
#ax.plot(xIAIM, f_exp(xIAIM, *poptIAIM), marker="", linestyle="-")

ax.legend()

plt.savefig("pt_time.pdf", bbox_inches="tight")

#plt.show()