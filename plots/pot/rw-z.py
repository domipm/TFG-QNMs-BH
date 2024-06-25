import numpy as np
import matplotlib.pyplot as plt

M = 1

def reggewheeler(x, l, s):
    
    return (1 - 2*M/x)*(l*(l+1)/x**2 + (1-s**2)*2*M/x**3)

def zerilli(x, l, s):

    return (1 - 2*M/x) * ( 2 * ((l-1)*(l+2)/2)**2 * ((l-1)*(l+2)/2 + 1) * x**3 + 6 * ( (l-1)*(l+2)/2 )**2 * M * x**2 + 18 * ((l-1)*(l+2)/2) * M**2 * x + 18 * M**3 ) / (x**3 * ( (l-1)*(l+2)/2 * x + 3*M )**2)


r = np.linspace(2.01,15,1000) # Schwarzschild coordinate

x = r + 2 * np.log(r/2 - 1)

plt.plot(x, reggewheeler(r, 2, 2), color="tab:blue", linestyle="-")
plt.plot(x, zerilli(r, 2, 2), color="tab:orange", linestyle="--")

plt.text(1,0.1,r"$\ell = 2$")
plt.text(1,0.32,r"$\ell = 3$")
plt.text(1,0.6,r"$\ell = 4$")

plt.plot(x, reggewheeler(r, 3, 2), color="tab:blue", linestyle="-", label=r"Regge-Wheeler $V_{\text{RW}}$")
plt.plot(x, zerilli(r, 3, 2), color="tab:orange", linestyle="--", label=r"Zerilli $V_{\text{Z}}$")

plt.plot(x, reggewheeler(r, 4, 2), color="tab:blue", linestyle="-")
plt.plot(x, zerilli(r, 4, 2), color="tab:orange", linestyle="--") # dashes=(5,4)

plt.title("Potenciales Regge-Wheeler y Zerilli")

plt.ylabel(r"Potencial $V(r)$")
plt.xlabel(r"Coordenada tortuga $r^*(r)$")

plt.legend()
plt.savefig("pot_comparison.pdf", bbox_inches="tight")
#plt.show()