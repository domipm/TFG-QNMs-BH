import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("conv.txt", skiprows=2)

plt.ylim(0.3, 0.4)

plt.hlines(0.3736716796, 2, 16, linestyle="--", color="black", zorder=20, label=r"Re($\omega_0$)")

plt.plot(data[:,0], data[:,1], label=r"-0.9")
plt.plot(data[:,0], data[:,3], label=r"-0.7")
plt.plot(data[:,0], data[:,5], label=r"-0.5")
plt.plot(data[:,0], data[:,7], label=r"-0.3")
plt.plot(data[:,0], data[:,9], label=r"-0.1")
#plt.plot(data[:,0], data[:,11], label=r"-0.01")
#plt.plot(data[:,0], data[:,13], label=r"-0.001")
plt.plot(data[:,0], data[:,15], label=r"0.000001")
#plt.plot(data[:,0], data[:,17], label=r"0.001")
plt.plot(data[:,0], data[:,19], label=r"0.1")
plt.plot(data[:,0], data[:,21], label=r"0.3")
plt.plot(data[:,0], data[:,23], label=r"0.5")
plt.plot(data[:,0], data[:,25], label=r"0.7")
plt.plot(data[:,0], data[:,27], label=r"0.9")

plt.legend(loc="upper left", bbox_to_anchor=(1, 0.85))
plt.savefig("conv_re.png", dpi=300, bbox_inches="tight")
#plt.show()

plt.close()

plt.ylim(-0.11, -0.05)

plt.hlines(-0.0889623141, 2, 16, linestyle="--", color="black", zorder=20, label=r"Im($\omega_0$)")

plt.plot(data[:,0], data[:,2], label=r"-0.9")
plt.plot(data[:,0], data[:,4], label=r"-0.7")
plt.plot(data[:,0], data[:,6], label=r"-0.5")
plt.plot(data[:,0], data[:,8], label=r"-0.3")
plt.plot(data[:,0], data[:,10], label=r"-0.1")
#plt.plot(data[:,0], data[:,11], label=r"-0.01")
#plt.plot(data[:,0], data[:,13], label=r"-0.001")
plt.plot(data[:,0], data[:,16], label=r"0.000001")
#plt.plot(data[:,0], data[:,17], label=r"0.001")
plt.plot(data[:,0], data[:,20], label=r"0.1")
plt.plot(data[:,0], data[:,22], label=r"0.3")
plt.plot(data[:,0], data[:,24], label=r"0.5")
plt.plot(data[:,0], data[:,26], label=r"0.7")
plt.plot(data[:,0], data[:,28], label=r"0.9")

plt.legend(loc="upper left", bbox_to_anchor=(1, 0.85))
plt.savefig("conv_im.png", dpi=300, bbox_inches="tight")
#plt.show()

plt.close()