import numpy as np
import matplotlib.pyplot as plt

M = 1

def reggewheeler(x, l, s):
    
    return (1 - 2*M/x)*(l*(l+1)/x**2 + (1-s**2)*2*M/x**3)

def zerilli(x, l, s):

    return (1 - 2*M/x) * ( 2 * ((l-1)*(l+2)/2)**2 * ((l-1)*(l+2)/2 + 1) * x**3 + 6 * ( (l-1)*(l+2)/2 )**2 * M * x**2 + 18 * ((l-1)*(l+2)/2) * M**2 * x + 18 * M**3 ) / (x**3 * ( (l-1)*(l+2)/2 * x + 3*M )**2)


r = np.linspace(2.01,15,1000) # Schwarzschild coordinate

x = r + 2 * np.log(r/2 - 1)

plt.plot(x, reggewheeler(r, 2, 0), color="tab:blue", linestyle="-.", label=r"$(s=0,l=2)$")
plt.plot(x, reggewheeler(r, 2, 1), color="tab:blue", linestyle="--", label=r"$(s=1,l=2)$")
plt.plot(x, reggewheeler(r, 2, 2), color="tab:blue", linestyle="-", label=r"$(s=2, l=2)$")

plt.title(r"Potenciales Regge-Wheeler perturbaciones spin $s$")

plt.ylabel(r"$V_{\text{eff}}$")
plt.xlabel(r"Coordenada $tortoise$ $r^*$")

plt.legend()
plt.savefig("pot_spin.pdf", bbox_inches="tight")
plt.show()