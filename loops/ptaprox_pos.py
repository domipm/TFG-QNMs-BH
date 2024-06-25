import sympy as sym
import numpy as np
import time
import sys

import aim_qnm_pos as aim

# Complex variable
I = sym.I
# Number of iterations to perform
n = 20

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

m = 1
l = 2
s = 2

# Fitted parameters -0.0925346788302107 0.14311474515072867
a = 0.16345200812517632
v = 0.15035326143963337

# Asymptotically Flat Schwarzschild Initial Parameters (lambda_0 and s_0)
l0 = 2*y*(1-I*w)/(1-y**2)
s0 = (1-2*I*w-2*w**2)/(2*(1-y**2))

# Find point where to substitute (maximum of potential)
#y0 = ( 1 - 2*m / ( 3*m/2 * 1/(l*(l+1)) * ( l*(l+1) - (1-s**2) + sym.sqrt( l**2*(l+1)**2 + 14/9*l*(l+1)*(1-s**2) + (1-s**2)**2 ) ) ) ).evalf(75)
y0 = 0

# AIM ALGORITHM

start = time.time()

# Create AIM solver object

aims = aim.aim_solver(l0,s0,x=y,x0=y0,n_iter=n)

a = np.linspace(-1, 1, 100)
print(a)

for p in a:

    print("\nPosition y0 = " + str(p) + "\n")

    aims.__init__(l0,s0,x=y,x0=p,n_iter=n)

    #aims.x0 = p

    aims.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
    aims.iaim_solve(solver="mpnum")

exit()

for pn in np.arange(-1.00,-0.01,0.01):

    print("\nPosition y0 = " + str(pn) + "\n")

    aims.x0 = pn

    aims.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
    aims.iaim_solve(solver="mpnum", print_delta = False)

for pp in np.arange(0.01,1.00,0.01):

    print("\nPosition y0 = " + str(pp) + "\n")

    aims.x0 = pp

    aims.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
    aims.iaim_solve(solver="mpnum", print_delta = False)