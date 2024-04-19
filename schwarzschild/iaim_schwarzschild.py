import sympy as sym
import numpy as np
import symengine as se
import mpmath as mp
import sys

# Set mpmath precision
mp.mp.dps = 50

# Complex variable
I = sym.I

# Symbolic variable definitions
y = sym.symbols("y", real=True)
w = sym.symbols("\omega")

# Angular momentum, spin, and mass
l = sym.symbols("l", real=True)
s = sym.symbols("s", real=True)
m = sym.symbols("m", real=True)

# Iterations to perform
N = 10

# Function to find maximum of potential (returns sympy float)
def max_point(m,l,s):

    return sym.S(1 - 2*m / ( 3*m/2 * 1/(l*(l+1)) * ( l*(l+1) - (1-s**2) + np.sqrt( l**2*(l+1)**2 + 14/9*l*(l+1)*(1-s**2) + (1-s**2)**2 ) ) ) )

print("Initializing parameters")

# Where to evaluate
y0 = max_point(1,2,0) 
# y0 = 0.322427767917328
y0 = 0 # Evaluating at 0 gets better results

print("y0 = " + str(y0) + "\n")

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

print("Computing Taylor series\n")

# Compute series expansion of l0 and s0 (works well, fast)

#l0_series = se.series(l0, y, y0, N).expand() # se.Add type

l0_series = se.series( se.series( l0, y, y0, N ), w, 0, N ).expand()

l0_coeff = np.zeros(N, dtype=object)
for i in range(0, N):
    #l0_coeff[i] = sym.nsimplify(sym.S((l0_series.expand().coeff(y,i)))) # Converts each coefficient to sympy and writes it as rational (exact?)
    l0_coeff[i] = l0_series.coeff(y,i).expand().evalf(N*N)
'''
eval_dict = {}
for a in l0_series.atoms(sym.Mul, sym.Add):
    numeric = [arg for arg in a.args if not arg.has(sym.Symbol)]
    symbolic = [arg for arg in a.args if arg.has(sym.Symbol)]
    eval_dict[a] = a.func(a.func(*numeric).evalf(N*N), a.func(*symbolic)) 
l0_series_eval = l0_series.subs(eval_dict)
'''
l0_series_eval = l0_series
print(l0_series_eval)
print("")

for i in range(0, N):
    l0_coeff[i] = l0_series_eval.coeff(y,i).expand()

print(l0_coeff)
print("")

#s0_series = se.series(s0, y, y0, N).expand() # se.Add type

s0_series = se.series( se.series( s0, y, y0, N ), w, 0, N ).expand()

s0_coeff = np.zeros(N, dtype=object)
for i in range(0, N):
    #s0_coeff[i] = sym.nsimplify(sym.S((s0_series.coeff(y,i)))) # Converts each coefficient to sympy and writes it as rational (exact?)
    s0_coeff[i] = s0_series.coeff(y,i).expand().evalf(N*N)
'''
eval_dict = {}
for a in s0_series.atoms(sym.Mul, sym.Add):
    numeric = [arg for arg in a.args if not arg.has(sym.Symbol)]
    symbolic = [arg for arg in a.args if arg.has(sym.Symbol)]
    eval_dict[a] = a.func(a.func(*numeric).evalf(N*N), a.func(*symbolic)) 
s0_series_eval = s0_series.subs(eval_dict)
'''
s0_series_eval = s0_series
for i in range(0, N):
    s0_coeff[i] = s0_series_eval.coeff(y,i).expand()

print(s0_series)
print("")

print(s0_coeff)

# Define matrices
C = np.zeros((N,N),dtype=object)
D = np.zeros((N,N),dtype=object)
# Initialize first column
C[:,0] = l0_coeff[:]
D[:,0] = s0_coeff[:]

#exit()

print("Calculating matrix coefficients\n") # Aquí ya hay algo que falla

'''
for n in range(0,N-1):
    for i in range(0,N-1-n):
        D[i,n+1] = ((i+1)*D[i+1,n]).expand()
        C[i,n+1] = ((i+1)*C[i+1,n]+D[i,n]).expand()
        for k in range(0,i+1):
            C[i,n+1] = (C[i,n+1] + C[k,0]*C[i-k,n]).expand()
            D[i,n+1] = (D[i,n+1] + D[k,0]*C[i-k,n]).expand()
'''

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

print("Computing characteristic polynomial\n")

# Compute polynomial d for last step
d = (D[0,N-1]*C[0,N-2] - D[0,N-2]*C[0,N-1]).expand()

print(d)

print( sym.nroots( sym.Poly(d,w) ) ) # Aquí ya hay algo que falla

d_coeff = sym.Poly(d.expand(),w).all_coeffs()

for i in range(len(d_coeff)):
    d_coeff[i] = mp.mpc(d_coeff[i])

print("Looking for solutions\n")

sols = mp.polyroots(d_coeff, maxsteps=100000, extraprec=15) # This one should be pretty exact (same one as for AIM method and it works perfectly)

for i in range(len(sols)): 
    sols[i] = mp.nstr(sols[i], 10)

print("Solutions found!\n")

f_sols = []
for i in range(len(sols)):
    if sym.re(sols[i]) > 0 and sym.im(sols[i]) < 0:
        f_sols.append(sols[i])

sols_sorted = sorted(f_sols, key = lambda x: sym.Abs(sym.im(x)))
for i in range(len(sols_sorted)):
    sym.pprint(sols_sorted[i])