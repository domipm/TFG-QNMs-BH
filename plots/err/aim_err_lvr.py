import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("err_leaver.txt", skiprows=1)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

plt.subplots_adjust(wspace=0, hspace=0)

axes[0].set_yscale("log")
axes[1].set_yscale("log")

# Parte real
axes[0].plot(data[:,0], data[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1.5)
axes[0].plot(data[:,0], data[:,3], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1.25)
axes[0].plot(data[:,0], data[:,5], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)
axes[0].plot(data[:,0], data[:,7], label=r"$\omega_3$", marker="x", linestyle="-", markersize=4, linewidth=0.75)
# Parte imaginaria
axes[1].plot(data[:,0], data[:,2], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1.5)
axes[1].plot(data[:,0], data[:,4], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1.25)
axes[1].plot(data[:,0], data[:,6], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)
axes[1].plot(data[:,0], data[:,8], label=r"$\omega_3$", marker="x", linestyle="-", markersize=4, linewidth=0.75)
#fig.tight_layout()

axes[0].set_ylim(10**-10, 10**0)
axes[1].set_ylim(10**-10, 10**0)

axes[0].set_title(r"Diferencia ref. Leaver p. real $\Delta{\text{Re}(\omega_n)}$")
axes[1].set_title(r"Diferencia ref. Leaver p. imaginaria $\Delta{\text{Im}(\omega_n)}$")

axes[1].yaxis.tick_right()

axes[0].set_ylabel(r"Diferencia $\log_{10}(\Delta\omega_{n})$")
axes[0].set_ylabel(r"Diferencia $\Delta\omega_{n}$")
#axes[1].set_ylabel(r"Error iterativo $\log_{10}(\varepsilon_{\text{Re}(\omega_n)})$")

axes[0].set_xlabel(r"Iteraciones $N$")
axes[1].set_xlabel(r"Iteraciones $N$")

axes[0].legend()
axes[1].legend()

plt.savefig("error_leaver.pdf", bbox_inches="tight")
plt.show()

exit()

error_re = np.loadtxt("aim_err.txt", skiprows=1)

plt.yscale("log")
plt.xlabel(r"Iteraciones")
plt.ylabel(r"Error Iterativo $\text{log}_{10}(\varepsilon_{\omega}$)")

plt.title(r"Error iterativo $\varepsilon_{\omega}$ de parte real Re($\omega_n$)")

plt.plot(error_re[:,0], error_re[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
plt.plot(error_re[:,0], error_re[:,2], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1)
plt.plot(error_re[:,0], error_re[:,3], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)
plt.tight_layout()
plt.legend()
plt.savefig("err_iter_re.pdf", dpi=300, bbox_inches='tight')
#plt.show()

plt.close()

#exit()

error_im = np.loadtxt("iter_error_im.txt", skiprows=1)

plt.yscale("log")
plt.xlabel(r"Iteraciones")
plt.ylabel(r"Error Iterativo $\varepsilon_{\omega}$")

plt.title(r"Error iterativo $\varepsilon_{\omega}$ de parte imaginaria Im($\omega_n$)")

plt.plot(error_im[:,0], error_im[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
plt.plot(error_im[:,0], error_im[:,2], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1)
plt.plot(error_im[:,0], error_im[:,3], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)
plt.tight_layout()
plt.legend()
plt.savefig("err_iter_im.pdf", dpi=300, bbox_inches='tight')
#plt.show()

plt.close()

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

plt.subplots_adjust(wspace=0, hspace=0)

axes[0].set_yscale("log")
axes[1].set_yscale("log")

axes[0].plot(error_re[:,0], error_re[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
axes[0].plot(error_re[:,0], error_re[:,2], label=r"$\omega_0$", marker="s", linestyle="-", markersize=4, linewidth=1)
axes[0].plot(error_re[:,0], error_re[:,3], label=r"$\omega_0$", marker="^", linestyle="-", markersize=4, linewidth=1)
axes[1].plot(error_im[:,0], error_im[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
axes[1].plot(error_im[:,0], error_im[:,2], label=r"$\omega_0$", marker="s", linestyle="-", markersize=4, linewidth=1)
axes[1].plot(error_im[:,0], error_im[:,3], label=r"$\omega_0$", marker="^", linestyle="-", markersize=4, linewidth=1)
#fig.tight_layout()

axes[0].set_ylim(10**-10, 10**2)
axes[1].set_ylim(10**-10, 10**2)

axes[0].set_title(r"Error iterativo parte real $\varepsilon_{\text{Re}(\omega_n)}$")
axes[1].set_title(r"Error iterativo parte imaginaria $\varepsilon_{\text{Im}(\omega_n)}$")

axes[1].yaxis.tick_right()

axes[0].set_ylabel(r"Error iterativo $\log_{10}(\varepsilon)$")
#axes[1].set_ylabel(r"Error iterativo $\log_{10}(\varepsilon_{\text{Re}(\omega_n)})$")

axes[0].set_xlabel("Iteraciones")
axes[1].set_xlabel("Iteraciones")

axes[0].legend()
axes[1].legend()

plt.savefig("error_iter.pdf", bbox_inches="tight")
plt.show()