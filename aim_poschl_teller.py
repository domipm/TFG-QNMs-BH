import numpy as np
import matplotlib.pyplot as plt

from sympy import *

#Number of iterations to perform
N = 5

#Where to evaluate the functions
y0 = 0

#Symbolic variable definitions
y = symbols("y", real=True)
w = symbols("\omega")

#Arrays for lambda and s parameters
l = np.empty(N,dtype=object)
s = np.empty(N,dtype=object)
#Arrays for lambda' and s' derivative of the parameters
lp = np.empty(N,dtype=object)
sp = np.empty(N,dtype=object)

#Initial parameters (specific to Poschl-Teller potential)
l[0] = (2*y*(1-I*w))/(1-y**2)
s[0] = (1-2*I*w-2*w**2)/(2*(1-y**2))

'''
AIM algorithm: 
-   Calculate derivative of previous parameters
    (lambda'_(n-1) and s'_(n-1))
-   Calculate new parameters
    (lambda_n and s_n in terms of previous
    params. and their derivatives)
-   Obtain delta (quantization condition)
-   Evaluate at point of interest (should be indep.)
-   Solve polynomial in w (omega) and filter sols.
'''

for n in range(1,N):

    print("\n*** ITERATION n=" + str(n) + " ***\n")

    #Calculate the previous derivatives
    lp[n-1] = diff(l[n-1],y)
    sp[n-1] = diff(s[n-1],y)

    #Calculate the new parameters lambda_n and s_n
    #Using the previous lambda_n-1, s_n-1, and derivatives
    #lambda'_n-1 and s'_n-1
    l[n] = lp[n-1] + s[n-1] + l[0]*l[n-1]
    s[n] = sp[n-1] + s[0]*l[n-1]

    #Quantization condition delta
    d = s[n]*l[n-1] - s[n-1]*l[n]

    #We find the "characteristic polynomial" by evaluating the 
    #Quantization condition and then simplify expression
    p = simplify(d.subs(y,y0))

    #We find all solutions to the polynomial in w (omega)
    sols = solve(p)
    print("All solutions: " + "\n" + str(sols) + "\n")
    #Filter solutions by taking the real part to be positive
    #And imaginary part to be negative (dampened oscillations)
    print("Filtered solutions:")
    for i in range(len(sols)):
        if re(sols[i]) > 0 and  im(sols[i]) < 0:
            print("w_" + str(i+1) + " = " + str(sols[i]))

    print("\n")

