import numpy as np
import time

import symengine as se
import mpmath as mp
import sympy as sym

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

#Precision of mpmath
mp.mp.dps = 50

start_T = time.time() #"Stopwatch" begin

#Number of iterations to perform 
#(add 2 with respect to AIM to obtain results up to same n index)
N = 6 # FAILS AT 8 ITERATIONS

#Symbolic variables definition
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

I = sym.I

#FUNCTION SERIES EXPANSION, REMOVE HIGHER ORDER, RETURN ARRAY WITH COEFFICIENTS
#Params: func. a, variable x, around point x0, order N (fixed)
def series_coeff(a,x,x0):

    a_series = se.series(a, x, x0, N).expand()
    coeff = np.array( sym.S( a_series.coeff(x,0).expand() ), dtype=object )
    for i in range(1,N):
        coeff = np.append(coeff, sym.S( a_series.coeff(x,i) ) )
    return coeff

#Initial parameters definition
l0 = (2*y*(1-I*w))/(1-y**2)
s0 = (1-2*I*w-2*w**2)/(2*(1-y**2))

y0 = 0

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

y0 = 0.322427767917328


#Matrix of coefficients needed
C = np.zeros((N,N),dtype=object)
D = np.zeros((N,N),dtype=object)
#First column initialized to coefficients of series expansion
C[:,0] = series_coeff(l0,y,y0)
D[:,0] = series_coeff(s0,y,y0)

#Calculate iteratively coefficients of C and D matrices
# PRECISION ERRORS MAIN SOURCE SHOULD BE HERE
# WE DRAG FLOATING POINT PRECISION ERRORS UP TO FINAL D POLYNOMIAL
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
d = (D[0,n]*C[0,n-1] - D[0,n-1]*C[0,n]).expand()

print(d)

#Find roots of the polynomial (algebraically) and display solutions/filtered solutions
#sols = sym.solve(d,w)                  #Algebraic root finder via sympy
#sols = sym.nroots(sym.Poly(d, w), n=8) #Numerical root finder via sympy

d_pol = sym.Poly(d, w)
dpol_coeff = d_pol.all_coeffs()
# Convert each coefficient into mpmath complex
for i in range(len(dpol_coeff)):
    dpol_coeff[i] = mp.mpc(dpol_coeff[i])
    
sols = mp.polyroots(dpol_coeff, maxsteps=10000, extraprec=150)

print("\nAll solutions:")
print(sols)
print("\nFiltered solutions:")
for i in range(len(sols)):
    if (sym.re(sols[i]) > 0 and  sym.im(sols[i]) < 0):
        print("w_" + str(i+1) + " = " + str(sols[i]))

end_T = time.time()#"Stopwatch" end
#Display total computation time
print("\nTotal computation time: " + str(end_T - start_T) + " s\n")