import sympy as sym
import time
import sys

sys.path.append('../')
import aim_qnm as aim

#Complex variable
I = sym.I
#Number of iterations to perform
n = 4
#Where to evaluate the functions
y0 = 0
#Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

#Poschl-Teller Initial Parameters (lambda_0 and s_0)
l0 = (2*y*(1-I*w))/(1-y**2)
s0 = (1-2*I*w-2*w**2)/(2*(1-y**2))

#AIM ALGORITHM

start = time.time()

aim = aim.aim_solver(l0,s0,y,y0,n)
aim.aim_init()
aim.aim_solve(solver="mpnum")

stop = time.time()

print("\n- Computation time (AIM Method): ", str(stop-start))

#IAIM ALGORITHM

start = time.time()

aim.iaim_init() #Initialize IAIM algorithm with same parameters as AIM class object
aim.iaim_solve(solver="mpnum", display_all=False, print_delta=False)

stop = time.time()

print("\n- Computation time (IAIM Method): ", str(stop-start) + "\n")