import sympy as sym
import numpy as np
import time
import sys

import aim_qnm_err as aim

# Complex variable
I = sym.I
# Number of iterations to perform
n = 64

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

# Angular momentum, spin, and mass
l = sym.symbols("l", real=True)
s = sym.symbols("s", real=True)
m = sym.symbols("m", real=True)

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

# Find point where to substitute (maximum of potential)
y0 = ( 1 - 2*m / ( 3*m/2 * 1/(l*(l+1)) * ( l*(l+1) - (1-s**2) + sym.sqrt( l**2*(l+1)**2 + 14/9*l*(l+1)*(1-s**2) + (1-s**2)**2 ) ) ) ).subs(m, m_val).subs(l, l_val).subs(s, s_val).evalf(75)

# AIM ALGORITHM

start = time.time()

# Create AIM solver object

aims = aim.aim_solver(l0,s0,x=y,x0=y0,n_iter=n)

aims.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
aims.iaim_solve(solver="mpnum")