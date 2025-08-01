import sys
import time
import numpy    as np
import sympy    as sym

# Import AIM library for computing QNMs
sys.path.append('../')
import aim_qnm  as aim

# Number of iterations to perform
n = 15

# Which potential to use (str : vec / sca / axi / pol)
pot = "vec"

# Numeric parameters (m=1, l=2)
m_val = 1
l_val = 2

# Parameter x_0 ( = 0 -> Schwarzschild )
x0_val = 10**-1
# Schwarzschild radius
rs_val = 1 * m_val

# Complex variable
I = sym.I

# Symbolic variable definitions
y, x = sym.symbols("y, x", real=True)
w = sym.symbols("\omega")

# Angular momentum, spin, and mass
l = sym.symbols("l", real=True)
s = sym.symbols("s", real=True)
m = sym.symbols("m", real=True)

# Schwarzschild radius
rs = sym.symbols("r_s", real = True)
# Parameter x_0
x0 = sym.symbols("x_0", real = True)


# Simpson-Visser parameters lambda_0 and s_0


# Scalar potential
if pot == "sca":
    # Simpson-Visser scalar parameters lambda_0 and s_0
    s0 = (16*m**4*w**2*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) + 16*m**4*w**2*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(-4*y**10 + 40*y**9 - 176*y**8 + 448*y**7 - 729*y**6 + 790*y**5 - 575*y**4 + 276*y**3 - 83*y**2 + 14*y - 1)*(1/(y**2 - 2*y + 1))**(13/2) + 8*I*m**3*w*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(2*y**9 - 18*y**8 + 69*y**7 - 149*y**6 + 201*y**5 - 177*y**4 + 103*y**3 - 39*y**2 + 9*y - 1)*(1/(y**2 - 2*y + 1))**(13/2) + 8*I*m**3*w*(y**2 - 2*y + 1)**6*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(-6*y**6 + 32*y**5 - 67*y**4 + 70*y**3 - 38*y**2 + 10*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) + 4*m**2*w**2*x0**2*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(4*y**10 - 40*y**9 + 176*y**8 - 448*y**7 + 729*y**6 - 790*y**5 + 575*y**4 - 276*y**3 + 83*y**2 - 14*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) + 4*m**2*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(-l**2*y**2 + 2*l**2*y - l**2 - l*y**2 + 2*l*y - l + y**4 - 4*y**3 + 8*y**2 - 8*y + 3)*(1/(y**2 - 2*y + 1))**(13/2) - 8*m**2*(y**2 - 2*y + 1)**6*(y**3 - 3*y**2 + 3*y - 1)*(-y**4 + 4*y**3 - 6*y**2 + 4*y - 1)**2*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) - 12*m**2*(y**2 - 2*y + 1)**2*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1) + 4*m**2*(y**2 - 2*y + 1)*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(l**2*y**2 - 2*l**2*y + l**2 + l*y**2 - 2*l*y + l - 2*y**2 + 4*y - 2) + 16*m**2*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)**2*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1) + 2*I*m*w*x0**2*(y**2 - 2*y + 1)**7*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(-2*y**9 + 18*y**8 - 69*y**7 + 149*y**6 - 201*y**5 + 177*y**4 - 103*y**3 + 39*y**2 - 9*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) + 2*I*m*w*x0**2*(y**2 - 2*y + 1)**6*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(8*y**8 - 58*y**7 + 178*y**6 - 301*y**5 + 305*y**4 - 188*y**3 + 68*y**2 - 13*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) + x0**2*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(-2*y**6 + 12*y**5 - 31*y**4 + 44*y**3 - 36*y**2 + 16*y - 3)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(1/(y**2 - 2*y + 1))**(13/2) + 4*x0**2*(y**2 - 2*y + 1)**3*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1) - x0**2*(y**2 - 2*y + 1)*(y**3 - 3*y**2 + 3*y - 1)*(-y**4 + 4*y**3 - 6*y**2 + 4*y - 1)**2*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1))/(y**2*(y - 1)**4*(-2*m + x0*y - x0)*(2*m + x0*y - x0)*(y**2 - 2*y + 1)**7*(y**3 - 3*y**2 + 3*y - 1)*(y**4 - 4*y**3 + 6*y**2 - 4*y + 1)*(y**5 - 5*y**4 + 10*y**3 - 10*y**2 + 5*y - 1)*(y**6 - 6*y**5 + 15*y**4 - 20*y**3 + 15*y**2 - 6*y + 1)*(1/(y**2 - 2*y + 1))**(13/2))
    l0 = -(32*I*m**3*w*y**2 - 64*I*m**3*w*y + 16*I*m**3*w - 12*m**2*y**2 + 16*m**2*y - 4*m**2 - 8*I*m*w*x0**2*y**4 + 32*I*m*w*x0**2*y**3 - 44*I*m*w*x0**2*y**2 + 24*I*m*w*x0**2*y - 4*I*m*w*x0**2 + 4*x0**2*y**4 - 13*x0**2*y**3 + 15*x0**2*y**2 - 7*x0**2*y + x0**2)/(y*(-4*m**2*y**2 + 8*m**2*y - 4*m**2 + x0**2*y**4 - 4*x0**2*y**3 + 6*x0**2*y**2 - 4*x0**2*y + x0**2))

    # Derivative of potential
    dv = x*(6*l**2*m*x**4 + 12*l**2*m*x**2*x0**2 + 6*l**2*m*x0**4 - 2*l**2*x**4*sym.sqrt(x**2 + x0**2) - 4*l**2*x**2*x0**2*sym.sqrt(x**2 + x0**2) - 2*l**2*x0**4*sym.sqrt(x**2 + x0**2) + 6*l*m*x**4 + 12*l*m*x**2*x0**2 + 6*l*m*x0**4 - 2*l*x**4*sym.sqrt(x**2 + x0**2) - 4*l*x**2*x0**2*sym.sqrt(x**2 + x0**2) - 2*l*x0**4*sym.sqrt(x**2 + x0**2) + 16*m**2*x**2*sym.sqrt(x**2 + x0**2) - 32*m**2*x0**2*sym.sqrt(x**2 + x0**2) - 6*m*x**4 + 18*m*x**2*x0**2 + 24*m*x0**4 - 4*x**2*x0**2*sym.sqrt(x**2 + x0**2) - 4*x0**4*sym.sqrt(x**2 + x0**2))/(x**2 + x0**2)**(9/2)
    # Search for maximum in x coordinate
    x_max = sym.nsolve(dv.subs( {rs: rs_val, l: l_val, x0: x0_val, m: m_val} ), x, 1)
    # Maximum point in y coordinate
    y_max = (1 - rs / sym.sqrt( x**2 + x0**2) ).subs( { x0: x0_val, rs: rs_val, x: x_max } )
    
    # Evaluation point
    y0 = y_max
    print("\ny_0 = ", y0, "\n")

# Vector potential
elif pot == "vec":
    # Simpson-Visser vector parameters lambda_0 and s_0
    s0 = 16*m**4*(-l**2*y**2/(4*m**2) + l**2*y**2/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) + l**2*y/(2*m**2) - l**2*y/(2*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - l**2/(4*m**2) + l**2/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - l*y**2/(4*m**2) + l*y**2/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) + l*y/(2*m**2) - l*y/(2*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - l/(4*m**2) + l/(4*m**2*sym.sqrt(1/(y**2 - 2*y + 1))) - 4*w**2*y**4 + 16*w**2*y**3 - 20*w**2*y**2 + 8*w**2*y - 2*I*w*y**4/m + 6*I*w*y**3/m - 6*I*w*y**2/m + 2*I*w*y/m + w**2*x0**2*y**6/m**2 - 6*w**2*x0**2*y**5/m**2 + 14*w**2*x0**2*y**4/m**2 - 16*w**2*x0**2*y**3/m**2 + 37*w**2*x0**2*y**2/(4*m**2) - 5*w**2*x0**2*y/(2*m**2) + w**2*x0**2/(4*m**2) + 3*I*w*x0**2*y**6/(4*m**3) - 15*I*w*x0**2*y**5/(4*m**3) + 59*I*w*x0**2*y**4/(8*m**3) - 57*I*w*x0**2*y**3/(8*m**3) + 27*I*w*x0**2*y**2/(8*m**3) - 5*I*w*x0**2*y/(8*m**3))/(y**2*(y - 1)**4*(-2*m + x0*y - x0)*(2*m + x0*y - x0))
    l0 = -(32*I*m**3*w*y**2 - 64*I*m**3*w*y + 16*I*m**3*w - 12*m**2*y**2 + 16*m**2*y - 4*m**2 - 8*I*m*w*x0**2*y**4 + 32*I*m*w*x0**2*y**3 - 44*I*m*w*x0**2*y**2 + 24*I*m*w*x0**2*y - 4*I*m*w*x0**2 + 4*x0**2*y**4 - 13*x0**2*y**3 + 15*x0**2*y**2 - 7*x0**2*y + x0**2)/(y*(-4*m**2*y**2 + 8*m**2*y - 4*m**2 + x0**2*y**4 - 4*x0**2*y**3 + 6*x0**2*y**2 - 4*x0**2*y + x0**2))

    # Evaluation point (constant)
    # y0 = 1.0 / 3.0
    y0 = 1 + 1 / x0_val

# Axial (tensor) potential
elif pot == "axi":
    # Simpson-Visser axial parameters lambda_0 and s_0
    s0 = (4*l**2*m**2*y**8*(1/(y**2 - 2*y + 1))**(5/2) - 32*l**2*m**2*y**7*(1/(y**2 - 2*y + 1))**(5/2) + 112*l**2*m**2*y**6*(1/(y**2 - 2*y + 1))**(5/2) - 224*l**2*m**2*y**5*(1/(y**2 - 2*y + 1))**(5/2) + 280*l**2*m**2*y**4*(1/(y**2 - 2*y + 1))**(5/2) - 224*l**2*m**2*y**3*(1/(y**2 - 2*y + 1))**(5/2) + 112*l**2*m**2*y**2*(1/(y**2 - 2*y + 1))**(5/2) - 4*l**2*m**2*y**2 - 32*l**2*m**2*y*(1/(y**2 - 2*y + 1))**(5/2) + 8*l**2*m**2*y + 4*l**2*m**2*(1/(y**2 - 2*y + 1))**(5/2) - 4*l**2*m**2 + 4*l*m**2*y**8*(1/(y**2 - 2*y + 1))**(5/2) - 32*l*m**2*y**7*(1/(y**2 - 2*y + 1))**(5/2) + 112*l*m**2*y**6*(1/(y**2 - 2*y + 1))**(5/2) - 224*l*m**2*y**5*(1/(y**2 - 2*y + 1))**(5/2) + 280*l*m**2*y**4*(1/(y**2 - 2*y + 1))**(5/2) - 224*l*m**2*y**3*(1/(y**2 - 2*y + 1))**(5/2) + 112*l*m**2*y**2*(1/(y**2 - 2*y + 1))**(5/2) - 4*l*m**2*y**2 - 32*l*m**2*y*(1/(y**2 - 2*y + 1))**(5/2) + 8*l*m**2*y + 4*l*m**2*(1/(y**2 - 2*y + 1))**(5/2) - 4*l*m**2 - 64*m**4*w**2*y**4 + 256*m**4*w**2*y**3 - 320*m**4*w**2*y**2 + 128*m**4*w**2*y - 32*I*m**3*w*y**4 + 96*I*m**3*w*y**3 - 96*I*m**3*w*y**2 + 32*I*m**3*w*y + 16*m**2*w**2*x0**2*y**6 - 96*m**2*w**2*x0**2*y**5 + 224*m**2*w**2*x0**2*y**4 - 256*m**2*w**2*x0**2*y**3 + 148*m**2*w**2*x0**2*y**2 - 40*m**2*w**2*x0**2*y + 4*m**2*w**2*x0**2 + 12*m**2*y**10*(1/(y**2 - 2*y + 1))**(7/2) - 120*m**2*y**9*(1/(y**2 - 2*y + 1))**(7/2) + 540*m**2*y**8*(1/(y**2 - 2*y + 1))**(7/2) - 1440*m**2*y**7*(1/(y**2 - 2*y + 1))**(7/2) + 2520*m**2*y**6*(1/(y**2 - 2*y + 1))**(7/2) - 3024*m**2*y**5*(1/(y**2 - 2*y + 1))**(7/2) + 2520*m**2*y**4*(1/(y**2 - 2*y + 1))**(7/2) - 12*m**2*y**4 - 1440*m**2*y**3*(1/(y**2 - 2*y + 1))**(7/2) + 48*m**2*y**3 + 540*m**2*y**2*(1/(y**2 - 2*y + 1))**(7/2) - 72*m**2*y**2 - 120*m**2*y*(1/(y**2 - 2*y + 1))**(7/2) + 48*m**2*y + 12*m**2*(1/(y**2 - 2*y + 1))**(7/2) - 12*m**2 + 12*I*m*w*x0**2*y**6 - 60*I*m*w*x0**2*y**5 + 118*I*m*w*x0**2*y**4 - 114*I*m*w*x0**2*y**3 + 54*I*m*w*x0**2*y**2 - 10*I*m*w*x0**2*y - 7*x0**2*y**10*(1/(y**2 - 2*y + 1))**(5/2) + 70*x0**2*y**9*(1/(y**2 - 2*y + 1))**(5/2) - 315*x0**2*y**8*(1/(y**2 - 2*y + 1))**(5/2) + 840*x0**2*y**7*(1/(y**2 - 2*y + 1))**(5/2) - 1470*x0**2*y**6*(1/(y**2 - 2*y + 1))**(5/2) + 4*x0**2*y**6 + 1764*x0**2*y**5*(1/(y**2 - 2*y + 1))**(5/2) - 24*x0**2*y**5 - 1470*x0**2*y**4*(1/(y**2 - 2*y + 1))**(5/2) + 63*x0**2*y**4 + 840*x0**2*y**3*(1/(y**2 - 2*y + 1))**(5/2) - 92*x0**2*y**3 - 315*x0**2*y**2*(1/(y**2 - 2*y + 1))**(5/2) + 78*x0**2*y**2 + 70*x0**2*y*(1/(y**2 - 2*y + 1))**(5/2) - 36*x0**2*y - 7*x0**2*(1/(y**2 - 2*y + 1))**(5/2) + 7*x0**2)/(y**2*(y - 1)**4*(-2*m + x0*y - x0)*(2*m + x0*y - x0))
    l0 = -(32*I*m**3*w*y**2 - 64*I*m**3*w*y + 16*I*m**3*w - 12*m**2*y**2 + 16*m**2*y - 4*m**2 - 8*I*m*w*x0**2*y**4 + 32*I*m*w*x0**2*y**3 - 44*I*m*w*x0**2*y**2 + 24*I*m*w*x0**2*y - 4*I*m*w*x0**2 + 4*x0**2*y**4 - 13*x0**2*y**3 + 15*x0**2*y**2 - 7*x0**2*y + x0**2)/(y*(-4*m**2*y**2 + 8*m**2*y - 4*m**2 + x0**2*y**4 - 4*x0**2*y**3 + 6*x0**2*y**2 - 4*x0**2*y + x0**2))

    # Derivative of potential
    dv = x*(6*l**2*m*x**2*sym.sqrt(x**2 + x0**2) + 6*l**2*m*x0**2*sym.sqrt(x**2 + x0**2) - 2*l**2*x**4 - 4*l**2*x**2*x0**2 - 2*l**2*x0**4 + 6*l*m*x**2*sym.sqrt(x**2 + x0**2) + 6*l*m*x0**2*sym.sqrt(x**2 + x0**2) - 2*l*x**4 - 4*l*x**2*x0**2 - 2*l*x0**4 - 48*m**2*x**2 + 48*m**2*x0**2 + 18*m*x**2*sym.sqrt(x**2 + x0**2) - 52*m*x0**2*sym.sqrt(x**2 + x0**2) + 12*x**2*x0**2 + 12*x0**4)/(x**2 + x0**2)**4
    # Search for maximum in x coordinate
    x_max = sym.nsolve(dv.subs( {rs: rs_val, l: l_val, x0: x0_val, m: m_val} ), x, 1)
    # Maximum point in y coordinate
    y_max = (1 - rs / sym.sqrt( x**2 + x0**2) ).subs( { x0: x0_val, rs: rs_val, x: x_max } )
    
    # Evaluation point
    y0 = y_max
    print("\ny_0 = ", y0, "\n")


# Polar (tensor) potential
elif pot == "pol":
    # Simpson-Visser polar parameters lambda_0 and s_0
    s0 = 0
    l0 = 0

    # Evaluation point (...)


# Substitute mass, ell, spin, and radius parameters
l0 = l0.subs(m,-m_val).subs(l,l_val).subs(x0, x0_val)
s0 = s0.subs(m,-m_val).subs(l,l_val).subs(x0, x0_val)


# AIM ALGORITHM


start = time.time()

# Create AIM solver object
aim = aim.aim_solver(l0, s0, x=y, x0=y0, n_iter=n)
# aim.aim_init() # Initialize AIM algorithm with same parameters as AIM class object
# aim.aim_solve(solver="mpnum", display_all=False, print_delta=False)

stop = time.time()

print("\nComputation time (AIM Method): ", str(stop-start) + "\n")


# IAIM ALGORITHM


start = time.time()

aim.iaim_init() # Initialize IAIM algorithm with same parameters as AIM class object
aim.iaim_solve(solver="mpnum")

stop = time.time()

print("\nComputation time (IAIM Method): ", str(stop-start) + "\n")