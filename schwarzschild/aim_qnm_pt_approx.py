import sympy as sym
import numpy as np
import time
import sys

sys.path.append('../')
import aim_qnm as aim

# Complex variable
I = sym.I
# Number of iterations to perform
n = 32

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

# Angular momentum, spin, and mass
l = sym.symbols("l", real=True)
s = sym.symbols("s", real=True)
m = sym.symbols("m", real=True)

# Numeric parameters (M=1, l=2, s=2)
m_val = 1
l_val = 2
s_val = 2

a = 1 / (m_val*sym.sqrt(27))
v = l_val*(l_val + 1)

# Asymptotically Flat Schwarzschild Initial Parameters (lambda_0 and s_0)
l0 =  2*y*(1-I*w)/(1-y**2)
s0 = (I*w-v)/(1-y**2) + (w**2*(1-a**2*y**2)/(a**2*(1-y**2)**2))

#Substitute
#l0 = l0.subs(m,m_val).subs(l,l_val).subs(s,s_val)
#s0 = s0.subs(m,m_val).subs(l,l_val).subs(s,s_val)

# Find point where to substitute (maximum of potential)
y0 = -( 1 - 2*m / ( 3*m/2 * 1/(l*(l+1)) * ( l*(l+1) - (1-s**2) + sym.sqrt( l**2*(l+1)**2 + 14/9*l*(l+1)*(1-s**2) + (1-s**2)**2 ) ) ) ).subs(m, m_val).subs(l, l_val).subs(s, s_val).evalf(75)

#y0 = 0.389

# Max value: 0.39038820320220757320583970795269124209880828857421875

# AIM ALGORITHM

start = time.time()

# Create AIM solver object
aim = aim.aim_solver(l0,s0,x=y,x0=y0,n_iter=n)
#aim.aim_init() # Initialize AIM algorithm with same parameters as AIM class object
#aim.aim_solve(solver="mpnum", display_all=False, print_delta=False)

stop = time.time()

print("\nComputation time (AIM Method): ", str(stop-start) + "\n")

# IAIM ALGORITHM (THIS ONE FAILS FOR NOW)

start = time.time()

aim.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
aim.iaim_solve(solver="mpnum", print_delta = False)

stop = time.time()

print("\nComputation time (IAIM Method): ", str(stop-start) + "\n")