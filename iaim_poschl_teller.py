from sympy import *
import numpy as np
import matplotlib.pyplot as plt

#NUMBER OF ITERATIONS TO PERFORM
N = 5

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

#print(series_coeff(l0,y,0))

#MATRIX OF COEFFICIENTS OF SERIES EXPANSION
C = np.zeros((N,N),dtype=object)
D = np.zeros((N,N),dtype=object)
#INITIALIZE FIRST COLUMN VALUES
C[:,0] = series_coeff(l0,y,0)
D[:,0] = series_coeff(s0,y,0)

for n in range(1,N):
    for i in range(N-1):
        C[i,n] = (i+1)*C[i+1,n-1] + D[i,n-1]
        D[i,n] = (i+1)*D[i+1,n-1]
        for k in range(0,i):
            C[i,n] = C[i,n] + C[k,0]*C[i-k,n-1]
            D[i,n] = D[i,n] + D[k,0]*C[i-k,n-1]

d = D[0,N-1]*C[0,N-2] - D[0,N-2]*C[0,N-1]

print(d)
print("\n")
print(solve(d,w))

exit()