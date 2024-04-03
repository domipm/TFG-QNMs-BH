import numpy as np
import time

from sympy import *

'''
IAIM (improved AIM) algorithm: 
-   Calculate series expansion of lambda_0 and s_0 coefficients
    around the point in which we will be evaluating
-   Create matrices C,D and store initial coefficients
-   Compute the necessary coefficients using the recursion relations
-   From the coefficients in the first row in the two last columns,
    apply the quantization condition
-   Solve the obtained polynomial in w (omega) algebraically
    and filter solutions
'''

start_T = time.time() #"Stopwatch" begin

#Number of iterations to perform 
#(add 2 with respect to AIM to obtain results up to same n index)
N = 6

#Symbolic variables definition
y = symbols("y", real=True)
w = symbols("\omega")

#FUNCTION SERIES EXPANSION, REMOVE HIGHER ORDER, RETURN ARRAY WITH COEFFICIENTS
#Params: func. a, variable x, around point x0, order N (fixed)
def series_coeff(a,x,x0):

    a_series = series(a,x,x0,N).removeO()
    coeff = np.array(a_series.subs(x,x0))
    for i in range(1,N):
        coeff = np.append(coeff, a_series.coeff(x**i))

    return coeff

#Initial parameters definition
l0 = (2*y*(1-I*w))/(1-y**2)
s0 = (1-2*I*w-2*w**2)/(2*(1-y**2))

#Matrix of coefficients needed
C = np.zeros((N,N),dtype=object)
D = np.zeros((N,N),dtype=object)
#First column initialized to coefficients of series expansion
C[:,0] = series_coeff(l0,y,0)
D[:,0] = series_coeff(s0,y,0)

#Calculate iteratively coefficients of C and D matrices
for n in range(0,N-1):
    for i in range(0,N):
        if (i+1 == N):
            D[i,n+1] = 0
            C[i,n+1] = D[i,n]
        else:
            D[i,n+1] = (i+1)*D[i+1,n]
            C[i,n+1] = (i+1)*C[i+1,n]+D[i,n]
        for k in range(0,i+1):
            C[i,n+1] = C[i,n+1] + C[k,0]*C[i-k,n]
            D[i,n+1] = D[i,n+1] + D[k,0]*C[i-k,n]

#Apply quantization condition to obtain polynomial in omega
d = D[0,n]*C[0,n-1] - D[0,n-1]*C[0,n]

#Find roots of the polynomial (algebraically) and display solutions/filtered solutions
sols = solve(d,w)
print("\nAll solutions:")
print(sols)
print("\nFiltered solutions:")
for i in range(len(sols)):
    if (re(sols[i]) > 0 and  im(sols[i]) < 0):
        print("w_" + str(i+1) + " = " + str(sols[i]))

end_T = time.time()#"Stopwatch" end
#Display total computation time
print("\nTotal computation time: " + str(end_T - start_T) + " s\n")