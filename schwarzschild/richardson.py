import numpy as np
import sympy as sym

'''
# IAIM RESULTS OBTAINED
w_0 = (0.37575090882991650 - 0.09024177409852527j)  n=2
w_0 = (0.37396196057878750 - 0.08888367557916661j) n=4
w_0 = (0.37366497335612225 - 0.08895537559014999j) n=8
w_0 = (0.37367171585027276 - 0.08896229704781046j) n=16
w_0 = (0.37367168439833615 - 0.08896231564583652j) n=32
'''

'''
# w_0
# Real part of quasinormal modes (Iterations n = 2, 4, 8, 16, 32, ...)
Rew = np.array([0.37575090882991650, 0.37396196057878750, 0.37366497335612225, 0.37367171585027276, 0.37367168439833615])
# Imaginary part of quasinormal modes (Iterations n = 2, 4, 8, 16, 32, ...)
Imw = np.array([-0.09024177409852527, -0.08888367557916661, -0.08895537559014999, -0.08896229704781046, -0.08896231564583652])
'''

'''
# w_2
Rew = np.array([0.22801169245848024, 0.27676301868174125, 0.24975442216089916])
Imw = np.array([-0.6841080802651178, -0.7148902235263268, -0.7039937637057415])
'''

# w_0 iter. (8,16,32,64)
# Using sym.Float does not change result
Rew = np.array([sym.Float('0.3739619605787875'), sym.Float('0.37366497335612225'), sym.Float('0.37367171585027276'), sym.Float('0.37367168439833615'), sym.Float('0.3736716846374484')])
Imw = np.array([-0.08888367557916661,-0.08895537559014999, -0.08896229704781046, -0.08896231564583652, -0.08896231576732479])

# Richardson extrapolation matrix (for real and imaginary parts)
D_re = np.zeros((len(Rew), len(Rew)))
D_im = np.zeros((len(Imw), len(Imw)))

D_re[:,0] = Rew
D_im[:,0] = Imw

# Interpolation for real and imaginary parts
for m in range(1, len(Rew)):
    for n in range(m, len(Rew)):
        D_re[n,m] = (4**m*D_re[n,m-1]-D_re[n-1,m-1])/(4**m-1)
        D_im[n,m] = (4**m*D_im[n,m-1]-D_im[n-1,m-1])/(4**m-1)

# Print lower-right term of each matrix
print("Re(w) = " + str( D_re[-1,-1] ))
print("Im(w) = " + str( D_im[-1,-1] ))

exit()

N = 3
d = np.zeros((N))
d = [Rew[0], Rew[1], Rew[2]]
D = np.zeros((N,N))
D[:,0] = d
for m in range(1,N):
    for n in range(m,N):
        D[n,m] = (4**m*D[n,m-1]-D[n-1,m-1])/(4**m-1)
D[2,2]