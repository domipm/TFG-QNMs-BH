#   Importing required libraries
import numpy as np      #   Used for numerical calculations
import sympy as sym     #   Used for symbolic calculations: Series expansion,
import time             #   Used to calculate computation times in critical parts
                        #   (better to do it in main code)
import symengine as se  #   Used for derivatives (faster than sympy)

import mpmath as mp     #   Used for precision handling

# Set mpmath precision and solver parameters
mp.mp.dps = 100

'''
Improved AIM algorithm (IAIM) taylored for computing quasinormal modes (QNMs)
of black holes (BHs).
Focused on solving second-order Schrödinger-like equation already written in the "AIM" form:

d^2/dx^2 psi(x) = lambda_0(x) d/dx psi(x) + s_0(x) psi(x)

Therefore, we have to compute the parameters lambda_0 and s_0 beforehand manually
(rewriting the potential in appropriate coordinates and including explicitly boundary conditions
into the solution, that is, incoming waves at horizon and outgoing at infinity)
'''

class aim_solver(object):

    w = sym.symbols("\omega") #  Symbol representing complex frequency of qnms

    #   Create AIM object with the necessary initial parameters, and create necessary arrays for lambda_n, s_n and derivatives
    #   Initial parameters depend only on symbols (x,y,...) defined in main code as necessary, with all numeric substitutions already performed
    def __init__(self, lambda_0, s_0, x, x0, n_iter): 

        #   Define initial parameters
        self.lambda_0 = lambda_0
        self.s_0 = s_0
        #   Number of iterations (size of arrays / matrices) to perform
        self.n_iter = n_iter
        #   Position coordinate and point of evaluation
        self.x = x
        self.x0 = x0
    
    '''
    AIM algorithm: 
    -   Begin by defining initial parameters lambda_0 and s_0
    -   Calculate necessary derivatives of "previous" parameters
        lambda'_(n-1) and s'_(n-1)
    -   Calculate new parameters lambda_n and s_n in terms of previous
        parameters and their derivatives
    -   Obtain delta (quantization condition)
    -   Evaluate at point of interest (maximum of potential)
    -   Solve polynomial in w (omega) and filter solutions
    '''
    
    #   Function to initialize all arrays needed for the AIM algorithm
    def aim_init(self):

        #   Arrays for lambda and s parameters
        self.l = np.empty(self.n_iter,dtype=object)
        self.s = np.empty(self.n_iter,dtype=object)
        #   Initialize array values
        self.l[0] = self.lambda_0
        self.s[0] = self.s_0
        #   Arrays for lambda' and s' derivative of the parameters
        self.lp = np.empty(self.n_iter,dtype=object)
        self.sp = np.empty(self.n_iter,dtype=object)
    
    #   Function to display results from AIM algorithm
    def aim_display(self, sols, display_all, n, n_modes = 100):

        print("\n*** AIM ITERATION n=" + str(n) + " ***\n")

        #   Display all solutions
        if (display_all == True): print("All solutions: " + "\n" + str(sols) + "\n")

        #   Filter solutions by positive real and negative imaginary part
        f_sols = []
        for i in range(len(sols)):
            if sym.re(sols[i]) > 0 and  sym.im(sols[i]) < 0:
                f_sols.append(sols[i])
        #   Display all filtered solutions
        if (display_all == True): print("Filtered solutions:\n" + str(f_sols) + "\n")

        #   Order solutions by imaginary part
        sols_sorted = sorted(f_sols, key = lambda x: sym.Abs(sym.im(x)))
        #   Display sorted solution by mode number (n)
        for i in range(len(sols_sorted)):
            if (i < n_modes): print("w_" + str(i) + " = " +  str(sols_sorted[i]))

    #   Solve via AIM algorithm
    def aim_solve(self, display_all = False, solver = "num", print_delta = False):

        w = sym.symbols("\omega") #  Symbol representing complex frequency of qnms

        for n in range(1, self.n_iter):

            #   Calculate the previous derivatives (via symengine)
            self.lp[n-1] = se.diff(self.l[n-1],self.x)
            self.sp[n-1] = se.diff(self.s[n-1],self.x)

            #   Calculate the new parameters lambda_n and s_n
            #   Using the previous lambda_n-1, s_n-1, and derivatives lambda'_n-1 and s'_n-1
            self.l[n] = self.lp[n-1] + self.s[n-1] + self.l[0]*self.l[n-1]
            self.s[n] = self.sp[n-1] + self.s[0]*self.l[n-1]

            #   Quantization condition delta / characteristic polynomial after substitution
            d = (self.s[n]*self.l[n-1] - self.s[n-1]*self.l[n]).subs(self.x,self.x0)
            #   "Normalize" quantization condition
            d = d / ( sym.Poly(d, self.w).all_coeffs() )[0]

            #   Display characteristic polynomial if needed
            if (print_delta == True): print(d.expand())

            #   Algebraic equation solver (via sympy)
            if (solver == "alg"): sols = sym.solve(d)
            #   Numeric polynomial root solver (via sympy)
            if (solver == "num"): sols = sym.nroots(sym.Poly(d, self.w), n=8, maxsteps=1000, cleanup=True)
            # Numerical polynomial root solver (via mpmath, faster and more precise)
            if (solver == "mpnum"):

                # Number of digits to evalf sympy expressions
                dps_evalf = 500
                # Maximum number of steps of solver
                nmax_solve = 10000
                # Extra precision of solver
                xprec_solve = 1000
                # Number of digits to print out
                dps_print = 10

                d_pol = sym.Poly(d, w)
                dpol_coeff = d_pol.all_coeffs()
                # Convert each coefficient into mpmath complex
                for i in range(len(dpol_coeff)):
                    # Evaluate sympy expression with 100 decimals gives precision up to 14 iterations
                    dpol_coeff[i] = mp.mpc(dpol_coeff[i].evalf(dps_evalf))

                sols = mp.polyroots(dpol_coeff, maxsteps=nmax_solve, extraprec=xprec_solve)

                # Write solutions with sig = 10 digits
                for i in range(len(sols)): 
                    sols[i] = mp.nstr(sols[i], dps_print)
                    
            #   Display the solution for each iteration
            self.aim_display(sols, display_all, n) 

    '''
    IAIM (improved AIM) algorithm: 
    -   Calculate series expansion of lambda_0 and s_0 coefficients
        around the point in which we will be evaluating
    -   Create matrices C,D and store initial coefficients
    -   Compute the necessary coefficients using the recursion relations
    -   From the coefficients in the first row in the two last columns,
        apply the quantization condition
    -   Solve the obtained polynomial in w (omega) and filter solutions
    ''' 

    #   FUNCTION SERIES EXPANSION, REMOVE HIGHER ORDER, RETURN ARRAY WITH COEFFICIENTS
    #   Params: func. a, variable x, around point x0, order N (fixed)
    def iaim_series_coeff(self,a,x,x0):

        print("Computing series expansion") # Debug
        start = time.time()

        #   Compute series via symengine
        a_series = se.series(a, x, x0, self.n_iter).expand()
        coeff = np.zeros(self.n_iter, dtype=object)
        for i in range(0,self.n_iter):
            coeff[i] = a_series.coeff(x,i)

        #   Compute series via sympy
        #a_series = sym.series(a,x,x0,self.n_iter).removeO().expand()
        #coeff = np.array(a_series.subs(x,x0).expand().evalf(500), dtype=object)
        #for i in range(1,self.n_iter):
        #    coeff = np.append(coeff, a_series.coeff(x**i).evalf(500))

        end = time.time()

        print(end - start)

        return coeff

    def iaim_init(self):

        #   Matrix of coefficients needed
        self.C = np.zeros((self.n_iter,self.n_iter),dtype=object)
        self.D = np.zeros((self.n_iter,self.n_iter),dtype=object)
        #   First column initialized to coefficients of series expansion
        self.C[:,0] = self.iaim_series_coeff(self.lambda_0,self.x,self.x0)
        self.D[:,0] = self.iaim_series_coeff(self.s_0,self.x,self.x0)
    
    def iaim_display(self, sols, display_all, n, n_modes = 100):

        print("\n*** IAIM ITERATION n=" + str(n+1) + " ***")

        #   Display all solutions
        if (display_all == True): print("\nAll solutions:" + str(sols))

        #   Filter solutions by positive real and negative imaginary part
        f_sols = []
        for i in range(len(sols)):
            if sym.re(sols[i]) > 0 and  sym.im(sols[i]) < 0:
                f_sols.append(sols[i])

        #   Display filtered (unordered) solutions
        if (display_all == True): print("\nFiltered solutions:" + str(f_sols) + "\n")

        #   Order solutions by imaginary part
        sols_sorted = sorted(f_sols, key = lambda x: sym.Abs(sym.im(x)))
        #   Display sorted solution by mode number (n)
        for i in range(len(sols_sorted)):
            if (i < n_modes): print("w_" + str(i) + " = " +  str(sols_sorted[i]))

    #   Solve via IAIM algorithm
    def iaim_solve(self, solver="num", display_all=False, print_delta=False):

        w = sym.symbols("\omega") #  Symbol representing complex frequency of qnms

        #   Compute iteratively coefficients / matrix elements
        for n in range(0,self.n_iter-1):

            # Debugging info
            print("Computing matrix elements for column " + str(n))
            
            for i in range(0,self.n_iter):
                if (i+1 == self.n_iter):
                    self.D[i,n+1] = 0
                    self.C[i,n+1] = (self.D[i,n]).expand()
                else:
                    self.D[i,n+1] = ((i+1)*self.D[i+1,n]).expand()
                    self.C[i,n+1] = ((i+1)*self.C[i+1,n]+self.D[i,n]).expand()
                for k in range(0,i+1):
                    self.C[i,n+1] = (self.C[i,n+1] + self.C[k,0]*self.C[i-k,n]).expand()
                    self.D[i,n+1] = (self.D[i,n+1] + self.D[k,0]*self.C[i-k,n]).expand()

        #   For each iteration compute the quantization condition and obtain modes
        for n in range(1, self.n_iter-1):

            #   Apply quantization condition and find characteristic polynomial
            d = (self.D[0,n]*self.C[0,n-1] - self.D[0,n-1]*self.C[0,n]).expand()

            if (print_delta == True): print(d.expand())

            #   Algebraic equation solver (via sympy)
            if (solver == "alg"): sols = sym.solve(d)

            #   Numeric polynomial root solver (via sympy)
            if (solver == "num"): sols = sym.nroots(sym.Poly(d, self.w),n=8, maxsteps=500, cleanup=True)

            # Numerical polynomial root solver (via mpmath, faster and more precise)
            if (solver == "mpnum"):

                # Number of digits to evalf sympy expressions
                dps_evalf = 100
                # Maximum number of steps of solver
                nmax_solve = 10000
                # Extra precision of solver
                xprec_solve = 100
                # Number of digits to print out
                dps_print = 10

                d_pol = sym.Poly(d, w)
                dpol_coeff = d_pol.all_coeffs()
                # Convert each coefficient into mpmath complex
                for i in range(len(dpol_coeff)):
                    # Evaluate sympy expression with 100 decimals gives precision up to 14 iterations
                    dpol_coeff[i] = mp.mpc(dpol_coeff[i].evalf(dps_evalf))

                sols = mp.polyroots(dpol_coeff, maxsteps=nmax_solve, extraprec=xprec_solve)

                # Write solutions with sig = 10 digits
                for i in range(len(sols)): 
                    sols[i] = mp.nstr(sols[i], dps_print)

            self.iaim_display(sols, display_all, n) #   Display the solution for each iteration