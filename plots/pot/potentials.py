import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

M = 1
L = 2
S = 2

v = L*(L+1) # Affect mainly height
a = 1/(M*np.sqrt(39.66)) # Affects height and spread

# sech^2(x) = 4 / ( e^-x + e^+x )^2

def poschlteller(x,x0):

    return 4*a**2*v / ( np.exp(-a*(x+x0)) + np.exp(a*(x+x0)) )**2

def poschlteller_fit(x,a0,v0):

    return 4*v0 / ( np.exp(-a0*(x-2.39)) + np.exp(a0*(x-2.39)) )**2

def reggewheeler(x):
    
    return (1 - 2*M/x)*(L*(L+1)/x**2 + (1-S**2)*2*M/x**3)

#plt.xlim(-1,10)

r = np.linspace(2.01,15,1000) # Schwarzschild coordinate

x = r + 2 * np.log(r/2 - 1)

vals = reggewheeler(r)

# We want to fit only values near maximum
vals_max = []
for i in range(len(vals)):
    if vals[i] > 0.13:
        vals_max = np.append(vals_max, vals[i])

indices = list(np.where(vals >= 0.13)[0])

print(len(vals_max))

popt, pcov = curve_fit(poschlteller_fit, x[indices], vals_max)

print(popt[0], popt[1])

#plt.title("Aproximación potencial Regge-Wheeler via Pöschl-Teller")
#plt.xlabel(r"Coordenada tortuga $r^*(r)$")
#plt.ylabel(r"Potencial $V(r)$")

plt.gca().set_aspect(200)

plt.plot(x, reggewheeler(r), label="Regge-Wheeler", linestyle="-", color="tab:blue")
plt.plot(x, poschlteller_fit(x, *popt), label="Pöschl-Teller", linestyle="--", color="tab:red")

#plt.plot(x, poschlteller(x, -2.39), label="Pöschl-Teller", linestyle="-")
#plt.legend()
plt.savefig("potential_approx_pres.pdf", bbox_inches="tight")
#plt.show()

exit()

print(reggewheeler(x))

rw_max = np.max(reggewheeler(x))
print(rw_max)

plt.plot(x_rw, reggewheeler(x_rw), label="R-W Pot")
plt.plot(x, poschlteller(x, -2.88), label="P-T Approx. Pot")
plt.legend()
plt.tight_layout()
plt.savefig("pot.png", dpi=300)
plt.show()