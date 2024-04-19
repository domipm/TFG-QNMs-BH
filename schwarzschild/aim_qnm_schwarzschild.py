import sympy as sym
import numpy as np
import time
import sys

sys.path.append('../')
import aim_qnm as aim

# Complex variable
I = sym.I
# Number of iterations to perform
n = 25

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

# Angular momentum, spin, and mass
l = sym.symbols("l", real=True)
s = sym.symbols("s", real=True)
m = sym.symbols("m", real=True)

# Function that returns the point of evaluation (xi_0 after change of variable)
def max_point(m,l,s):

    return sym.S(1 - 2*m / ( 3*m/2 * 1/(l*(l+1)) * ( l*(l+1) - (1-s**2) + np.sqrt( l**2*(l+1)**2 + 14/9*l*(l+1)*(1-s**2) + (1-s**2)**2 ) ) ) )

# Where to evaluate
y0 = max_point(1,2,0)

# Asymptotically Flat Schwarzschild Initial Parameters (lambda_0 and s_0)
l0 = (4*m*I*w*(2*y**2 - 4*y + 1) - (1 - 3*y)*(1-y) ) / (y*(1-y)**2) 
s0 = ( 16*m**2*w**2*(y-2) - 8*m*I*w*(1-y) + l*(l+1) + (1-s**2)*(1-y) ) / ( y*(1-y)**2 )

# Numeric parameters (M=1, l=2, s=2)
m_val = 1
l_val = 2
s_val = 2
#Substitute
l0 = l0.subs(m,m_val).subs(l,l_val).subs(s,s_val)
s0 = s0.subs(m,m_val).subs(l,l_val).subs(s,s_val)

# AIM ALGORITHM

start = time.time()

# Create AIM solver object
aim = aim.aim_solver(l0,s0,x=y,x0=y0,n_iter=n)
aim.aim_init() # Initialize AIM algorithm with same parameters as AIM class object
#aim.aim_solve(solver="mpnum", display_all=False, print_delta=False)

stop = time.time()

print("\nComputation time (AIM Method): ", str(stop-start) + "\n")

# IAIM ALGORITHM (THIS ONE FAILS FOR NOW)

start = time.time()

aim.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
aim.iaim_solve(solver="mpnum", print_delta = False)

stop = time.time()

print("\nComputation time (IAIM Method): ", str(stop-start) + "\n")