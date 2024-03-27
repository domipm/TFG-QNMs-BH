import numpy as np
import matplotlib.pyplot as plt

from sympy import *

#Number of iterations to perform
N = 10

#Symbolic variable definitions
y = symbols("y", real=True)
w = symbols("\omega")

l = np.empty(N,dtype=object)
s = np.empty(N,dtype=object)

lp = np.empty(N,dtype=object)
sp = np.empty(N,dtype=object)

#Initial parameters (problem specific)
l[0] = (2*y*(1-I*w))/(1-y**2)
s[0] = (1-2*I*w-2*w**2)/(2*(1-y**2))

lp[0] = diff(l[0],y)
sp[0] = diff(s[0],y)

#AIM algorithm, calculate derivative of previous
#parameters, calculate new parameters
#obtain delta, evaluate, solve

for n in range(1,N):

    print("\n*** ITERATION " + str(n) + " ***\n")

    lp[n] = diff(l[n-1],y)
    sp[n] = diff(s[n-1],y)

    l[n] = lp[n-1] + s[n-1] + l[0]*l[n-1]
    s[n] = sp[n-1] + s[0]*l[n-1]

    d = s[n]*l[n-1] - s[n-1]*l[n]

    p = simplify(d.subs(y,0))

    sols = solve(p)
    print("All solutions: " + "\n" + str(sols) + "\n")
    print("Filtered solutions:")
    for i in range(len(sols)):
        if re(sols[i]) > 0 and  im(sols[i]) < 0:
            print("w_" + str(i+1) + " = " + str(sols[i]))

    print("\n")

