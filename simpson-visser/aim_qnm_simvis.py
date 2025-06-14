import sympy as sym
import numpy as np
import time
import sys

sys.path.append('../')
import aim_qnm as aim

# Complex variable
I = sym.I
# Number of iterations to perform
n = 5

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

# Angular momentum, spin, and mass
l = sym.symbols("l", real=True)
s = sym.symbols("s", real=True)
m = sym.symbols("m", real=True)

# Schwarzschild radius
rs = sym.symbols("r_s", real = True)
# Parameter x_0
x0 = sym.symbols("x_0", real = True)

# Simpson-Visser parameters lambda_0 and s_0
s0 = 16*m**4*(-l**2*y**2/(4*m**2) + l**2*y**2/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) + l**2*y/(2*m**2) - l**2*y/(2*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - l**2/(4*m**2) + l**2/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - l*y**2/(4*m**2) + l*y**2/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) + l*y/(2*m**2) - l*y/(2*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - l/(4*m**2) + l/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - 4*w**2*y**4 + 16*w**2*y**3 - 20*w**2*y**2 + 8*w**2*y - 2*I*w*y**4/m + 6*I*w*y**3/m - 6*I*w*y**2/m + 2*I*w*y/m + w**2*x0**2*y**6/m**2 - 6*w**2*x0**2*y**5/m**2 + 14*w**2*x0**2*y**4/m**2 - 16*w**2*x0**2*y**3/m**2 + 37*w**2*x0**2*y**2/(4*m**2) - 5*w**2*x0**2*y/(2*m**2) + w**2*x0**2/(4*m**2) + 3*I*w*x0**2*y**6/(4*m**3) - 15*I*w*x0**2*y**5/(4*m**3) + 59*I*w*x0**2*y**4/(8*m**3) - 57*I*w*x0**2*y**3/(8*m**3) + 27*I*w*x0**2*y**2/(8*m**3) - 5*I*w*x0**2*y/(8*m**3))/(y**2*(y - 1)**4*(-2*m + x0*y - x0)*(2*m + x0*y - x0))
l0 = (-32*I*m**3*w*y**4 + 128*I*m**3*w*y**3 - 176*I*m**3*w*y**2 + 96*I*m**3*w*y - 16*I*m**3*w + 12*m**2*y**4 - 40*m**2*y**3 + 48*m**2*y**2 - 24*m**2*y + 4*m**2 + 8*I*m*w*x0**2*y**6 - 48*I*m*w*x0**2*y**5 + 116*I*m*w*x0**2*y**4 - 144*I*m*w*x0**2*y**3 + 96*I*m*w*x0**2*y**2 - 32*I*m*w*x0**2*y + 4*I*m*w*x0**2 - 4*x0**2*y**6 + 21*x0**2*y**5 - 45*x0**2*y**4 + 50*x0**2*y**3 - 30*x0**2*y**2 + 9*x0**2*y - x0**2)/(y*(y - 1)**4*(-2*m + x0*y - x0)*(2*m + x0*y - x0))

# Numeric parameters (m=1, l=2)
m_val = 1
l_val = 2

# Parameter x_0 ( = 0 -> Schwarzschild )
x0_val = 0.1

# Substitute mass, ell, spin, and radius parameters
l0 = l0.subs(m,m_val).subs(l,l_val).subs(x0, x0_val)
s0 = s0.subs(m,m_val).subs(l,l_val).subs(x0, x0_val)

# Find point where to substitute (maximum of potential)
y0 = 1/3

# AIM ALGORITHM

start = time.time()

# Create AIM solver object
aim = aim.aim_solver(l0, s0, x=y, x0=y0, n_iter=n)
aim.aim_init() # Initialize AIM algorithm with same parameters as AIM class object
aim.aim_solve(solver="mpnum", display_all=False, print_delta=False)

stop = time.time()

print("\nComputation time (AIM Method): ", str(stop-start) + "\n")

# IAIM ALGORITHM

start = time.time()

aim.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
aim.iaim_solve(solver="mpnum", print_delta = False)

stop = time.time()

print("\nComputation time (IAIM Method): ", str(stop-start) + "\n")