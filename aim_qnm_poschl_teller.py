import aim_qnm as aim
import sympy as sym
import time

#Complex variable
I = sym.I
#Number of iterations to perform
n = 5
#Where to evaluate the functions
y0 = 0
#Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

#Poschl-Teller Initial Parameters (lambda_0 and s_0)
l0 = (2*y*(1-I*w))/(1-y**2)
s0 = (1-2*I*w-2*w**2)/(2*(1-y**2))

start = time.time()

aim = aim.aim_solver(l0,s0,n)
aim.aim_init()
aim.aim_solve(display_all=True,solver="alg",x=y,x0=y0)

stop = time.time()

print("\nComputation time (AIM Method): ", str(stop-start) + "\n")