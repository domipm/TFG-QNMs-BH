import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.txt", skiprows=1)

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

plt.plot(data[:,0], data[:,2], label=r"$\omega_0$", marker="o", linestyle="-", markersize=4, linewidth=1)
plt.plot(data[:,0], data[:,4], label=r"$\omega_1$", marker="s", linestyle="-", markersize=4, linewidth=1)
plt.plot(data[:,0], data[:,6], label=r"$\omega_2$", marker="^", linestyle="-", markersize=4, linewidth=1)

plt.vlines(x=0.3903882032022075, ymin=-0.75, ymax=1.25, linestyle="--", color="grey", label=r"$y_{\text{max}}$")
plt.vlines(x=0, ymin=-0.75, ymax=1.25, linestyle="--", color="grey", label=r"$y_{\text{EH}}$")

plt.xlabel(r"y")
plt.ylabel(r"Im($\omega_n$)")

plt.tight_layout()
plt.legend()
plt.savefig("pos_one_im.png", dpi=300, bbox_inches="tight")
#plt.show()

plt.close()

f, (ax_re, ax_im) = plt.subplots(2, 1, sharex=True, facecolor='w')

plt.subplots_adjust(wspace=0, hspace=0)

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

ax_re.vlines(x=0.3903882032022075, ymin=-1.55, ymax=1.25, linestyle="--", color="dimgrey", label=r"$y_{\text{max}}$", linewidth=0.95)
ax_re.vlines(x=0, ymin=-2.55, ymax=1.75, linestyle="-", color="dimgrey", label=r"$y_{\text{EH}}$", linewidth=0.95)
ax_im.vlines(x=0.3903882032022075, ymin=-1.55, ymax=1.25, linestyle="--", color="dimgrey", label=r"$y_{\text{max}}$", linewidth=0.95)
ax_im.vlines(x=0, ymin=-2.55, ymax=1.75, linestyle="-", color="dimgrey", label=r"$y_{\text{EH}}$", linewidth=0.95)

ax_re.set_ylim(-0.25,1.25)
ax_im.set_ylim(-0.75,1.25)

#ax_re.spines['bottom'].set_visible(False)
#ax_im.spines['top'].set_visible(False)
#ax_re.xaxis.tick_bottom()
#ax_re.tick_params(labeltop='off')  # don't put tick labels at the top
#ax_im.xaxis.tick_bottom()

ax_re.legend(loc = "upper left", bbox_to_anchor=(1, 0.32))
#ax_im.legend(loc = "upper left")

ax_re.set_title("Frecuencias en función del punto de evaluación")
ax_im.set_xlabel(r"$y_0$")
#ax_re.set_xlabel(r"$y_0$")

plt.savefig("pos_one_split.pdf", dpi=300, bbox_inches="tight")
plt.show()