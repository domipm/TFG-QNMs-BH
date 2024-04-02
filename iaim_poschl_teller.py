from sympy import *
import numpy as np
import matplotlib.pyplot as plt

#NUMBER OF ITERATIONS TO PERFORM
N = 15

#SYMBOLIC VARIABLES DEFINITIONS
y = symbols("y", real=True)
w = symbols("\omega")

#FUNCTION SERIES EXPANSION, REMOVE O, RETURN ARRAY WITH COEFFICIENTS
#Params: func. a, variable x, around point x0, order N (fixed)
def series_coeff(a,x,x0):

    a_series = series(a,x,x0,N).removeO()
    coeff = np.array(a_series.subs(x,x0))
    for i in range(1,N):
        coeff = np.append(coeff, a_series.coeff(x**i))

    return coeff

#INITIAL PARAMETERS
l0 = (2*y*(1-I*w))/(1-y**2)
s0 = (1-2*I*w-2*w**2)/(2*(1-y**2))

#MATRIX OF COEFFICIENTS OF SERIES EXPANSION
C = np.zeros((N,N),dtype=object)
D = np.zeros((N,N),dtype=object)
#INITIALIZE FIRST COLUMN VALUES
C[:,0] = series_coeff(l0,y,0)
D[:,0] = series_coeff(s0,y,0)

#CALCULATE MATRICES C and D
for n in range(0,N-1): #COLUMNS
    for i in range(0,N): #ROWS
        if (i+1 == N):
            D[i,n+1] = 0
            C[i,n+1] = D[i,n]
        else:
            D[i,n+1] = (i+1)*D[i+1,n]
            C[i,n+1] = (i+1)*C[i+1,n]+D[i,n]
        for k in range(0,i+1): #SUM INTERIOR
            C[i,n+1] = C[i,n+1] + C[k,0]*C[i-k,n]
            D[i,n+1] = D[i,n+1] + D[k,0]*C[i-k,n]

#SOLVE QUANTIZATION CONDITION
d = D[0,n]*C[0,n-1] - D[0,n-1]*C[0,n]
#d.expand()

#FIND ROOTS OF POLYNOMIAL OF w (omega) (Using sympy, algebraic roots, might be slow, use numerical instead)
sols = solve(d,w)
print("All solutions:")
print(sols)
print("Filtered solutions:")
for i in range(len(sols)):
    if (re(sols[i]) > 0 and  im(sols[i]) < 0):
        print("w_" + str(i+1) + " = " + str(sols[i]))

exit()