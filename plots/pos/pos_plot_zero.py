import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("near_zero.txt", skiprows=1)

plt.plot(data[:,0], data[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
plt.plot(data[:,0], data[:,3], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1)
plt.plot(data[:,0], data[:,5], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)

plt.vlines(x=0.3903882032022075, ymin=-0.25, ymax=1.25, linestyle="--", color="grey", label=r"$y_{\text{max}}$")
plt.vlines(x=0, ymin=-0.25, ymax=1.25, linestyle="--", color="grey", label=r"$y_{\text{EH}}$")

plt.xlabel(r"y")
plt.ylabel(r"Re($\omega_n$)")

plt.tight_layout()
plt.legend()
plt.savefig("pos_one_re.png", dpi=300, bbox_inches="tight")
#plt.show()

plt.close()

f, (ax_re, ax_im) = plt.subplots(2, 1, sharex=True, facecolor='w')

#ax_re.set_title("real")
#ax_im.set_title("imag")

#ax_re.get_xaxis().set_visible(False)

ax_re.set_ylabel(r"Re($\omega_n$)")
ax_im.set_ylabel(r"Im($\omega_n$)")

ax_re.plot(data[:,0], data[:,1], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
ax_re.plot(data[:,0], data[:,3], label=r"$\omega_0$", marker="s", linestyle="-", markersize=4, linewidth=1)
ax_re.plot(data[:,0], data[:,5], label=r"$\omega_0$", marker="^", linestyle="-", markersize=4, linewidth=1)

ax_im.plot(data[:,0], data[:,2], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
ax_im.plot(data[:,0], data[:,4], label=r"$\omega_0$", marker="s", linestyle="-", markersize=4, linewidth=1)
ax_im.plot(data[:,0], data[:,6], label=r"$\omega_0$", marker="^", linestyle="-", markersize=4, linewidth=1)

#ax_re.vlines(x=0.3903882032022075, ymin=-1.55, ymax=1.25, linestyle="--", color="dimgrey", label=r"$y_{\text{max}}$", linewidth=0.95)
#ax_re.vlines(x=0, ymin=-2.55, ymax=1.75, linestyle="-", color="dimgrey", label=r"$y_{\text{EH}}$", linewidth=0.95)
#ax_im.vlines(x=0.3903882032022075, ymin=-1.55, ymax=1.25, linestyle="--", color="dimgrey", label=r"$y_{\text{max}}$", linewidth=0.95)
#ax_im.vlines(x=0, ymin=-2.55, ymax=1.75, linestyle="-", color="dimgrey", label=r"$y_{\text{EH}}$", linewidth=0.95)

ax_re.set_ylim(-0.25,3.25)
ax_im.set_ylim(-1.25,1.25)

ax_re.set_xlim(10e-7, 10e-2)

ax_re.set_xscale("log")
ax_im.set_xscale("log")

#ax_re.spines['bottom'].set_visible(False)
#ax_im.spines['top'].set_visible(False)
#ax_re.xaxis.tick_bottom()
#ax_re.tick_params(labeltop='off')  # don't put tick labels at the top
#ax_im.xaxis.tick_bottom()

ax_re.legend(loc = "upper left", bbox_to_anchor=(1, 0.12))
#ax_im.legend(loc = "upper left")

ax_re.set_title("Frecuencias en función del punto de evaluación")
ax_im.set_xlabel(r"$y_0$")
#ax_re.set_xlabel(r"$y_0$")

plt.savefig("pos_zero.png", dpi=300, bbox_inches="tight")
#plt.show()

exit()

ax_re.set_title(r"Estructura de bandas nanowire ")
ax_im.set_xlabel(r"$k$ [nm$^{-1}$]")
ax_re.get_xaxis().set_visible(False)
ax_re.set_ylabel(r"$E$ [eV]")
ax_re.yaxis.set_label_coords(0.05, 0.5, transform=f.transFigure)
# Máximos y mínimos de bandas
ax_im.hlines(max, k[0,0], k[0,-1], linestyle="--", alpha=0.65, color="black", label="BV max = " + str(round(max,5)) + " [eV]")
ax_re.hlines(min, k[0,0], k[0,-1], linestyle="--", alpha=0.65, color="black", label="BC min = " + str(round(min,5)) + " [eV]")
# Todas las bandas
for i in range(0, e_prec):
    ax_re.plot(k[i], e[i], marker=".", markersize=3, linestyle="-")
    ax_im.plot(k[i], e[i], marker=".", markersize=3, linestyle="-")
# Find limits of Y values
ymax = np.max(e[:,:])
ymin = np.min(e[:,:])
# Limites eje-Y
ax_re.set_ylim(min-0.05, np.max(e[:,:])+0.05)  # Banda de conducción
ax_im.set_ylim(np.min(e[:,:])-0.05, max+0.05)  # Banda de valencia
# hide the spines between ax and ax2
ax_re.spines['bottom'].set_visible(False)
ax_im.spines['top'].set_visible(False)
ax_re.xaxis.tick_top()
ax_re.tick_params(labeltop='off')  # don't put tick labels at the top
ax_im.xaxis.tick_bottom()
d = .005  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax_re.transAxes, color='k', clip_on=False)
ax_re.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax_re.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax_im.transAxes)  # switch to the bottom axes
ax_im.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax_im.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
plt.savefig("band_structure_split_"  + str(crystal) + "_" + str(D) + ".png", dpi=300, bbox_inches="tight")
#plt.show()

    

exit()

plt.plot(data[:,0], data[:,2], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
plt.plot(data[:,0], data[:,4], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1)
plt.plot(data[:,0], data[:,6], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)

plt.vlines(x=0.3903882032022075, ymin=-0.25, ymax=1.25, linestyle="--", color="grey", label=r"$y_{\text{max}}$")
plt.vlines(x=0, ymin=-0.25, ymax=1.25, linestyle="--", color="grey", label=r"$y_{\text{EH}}$")

plt.xlabel(r"y")
plt.ylabel(r"Im($\omega_n$)")

plt.tight_layout()
plt.legend()
plt.savefig("pos_one_im.png", dpi=300, bbox_inches="tight")
plt.show()

plt.close()