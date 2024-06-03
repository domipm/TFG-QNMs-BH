import sympy as sym
import numpy as np
import time
import sys

sys.path.append('../')
import aim_qnm as aim

# Complex variable
I = sym.I
# Number of iterations to perform
n = 20

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

# Fitted parameters
# ( Correct parameters checked in ref. )
a = 0.16345200812517632
v = 0.15035326143963337

#a = 1
#v = 1/2

# Asymptotically Flat Schwarzschild Initial Parameters (lambda_0 and s_0)
l0 = (2*y*(1-I*w))/(1-y**2)
s0 = (v-w*(w+I*a))/(a**2*(1-y**2))

# Find point where to substitute
y0 = 0

# AIM ALGORITHM

start = time.time()

# Create AIM solver object ( Gives same results as IAIM )
aim = aim.aim_solver(l0,s0,x=y,x0=y0,n_iter=n)
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