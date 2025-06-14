import  numpy       as np       #   Used for numerical calculations
import  sympy       as sym      #   Used for symbolic calculations: series expansion
import  mpmath      as mp       #   Used for precision handling
import  symengine   as se       #   Used for derivatives (faster than sympy)


'''

Implementation of the (Improved) AIM algorithm taylored for computing quasinormal modes (QNMs)
of black holes (BHs) geometries. Used for solving the second-order Schrödinger-like equation of the form:

d^2/dx^2 psi(x) = lambda_0(x) d/dx psi(x) + s_0(x) psi(x)

Requires previous computation of the parameters lambda_0 and s_0 in terms of coordinate x, including implicitly
the boundary conidtions, to be evaluated for n_iter at the evaluation point x_0.

'''


#   Symbol representing complex frequency of qnms
w = sym.symbols("\omega")


class aim_solver(object):


    #   Create AIM object with the necessary initial parameters, and create necessary arrays for lambda_n, s_n and derivatives
    #   Initial parameters depend only on symbols (x,y,...) defined in main code as necessary, with all numeric substitutions already performed
    def __init__(self, lambda_0, s_0, x, x0, n_iter): 

        #   Define initial parameters
        self.lambda_0 = lambda_0
        self.s_0 = s_0
        #   Number of iterations (size of arrays / matrices) to perform
        self.n_iter = n_iter+1
        #   Position coordinate and point of evaluation
        self.x = x
        self.x0 = x0


    #   Function to display results from (I)AIM algorithm
    def aim_display(self, which, sols, display_all, n, n_modes = 15):
        '''
        Auxiliary function to display results on screen. If no stable solutions are found (with Im < 0),
        unstable solutions are included (with Im > 0).

        Args:
            which (str): Which algorithm is being used AIM/IAIM
            sols (list): Found solutions to the polynomial equation in omega
            display_all (bool): Whether to show all solutions (positive and negative imaginary part)
            n (int): Iteration number
            n_modes (int): Maximum number of frequencies to display
        
        Returns:
            sols_sorted: Sorted solutions / quasinormal frequencies
        '''

        print(f"\n*** {which} ITERATION n={str(n)} ***\n")

        # Display only filtered solutions
        if (display_all == False):

            #   Filter solutions by positive real and negative imaginary part
            f_sols = []
            # Go over each solution
            for i in range(len(sols)):
                # Check if real part is positive and imaginary part is negative
                if sym.re(sols[i]) > 0 and sym.im(sols[i]) < 0:
                    f_sols.append(sols[i])

            # Check if no stable solutions are found
            if len(f_sols) == 0:
                # Include unstable solutions
                for i in range(len(sols)):
                    # Check if both real and imaginary part are positive
                    if sym.re(sols[i]) > 0 and sym.im(sols[i]) > 0:
                        # Append unstable solution
                        f_sols.append(sols[i])

        # Otherwise, display all solutions (no filtering)
        elif (display_all == True):
            # Consider all solutions
            f_sols = sols

        #   Order solutions by imaginary part
        sols_sorted = sorted(f_sols, key = lambda x: sym.Abs(sym.im(x)))
        #   Display sorted solution by mode number (n)
        for i in range(len(sols_sorted)):
            if (i < n_modes): 
                print( str( float(mp.re(sols_sorted[i] ))) + "\t" + str( float(mp.im(sols_sorted[i]) )) )  

        # Return filtered solutions
        return sols_sorted
    

    #   Function to solve polynomial 
    def aim_poly_solve(self, d, solver, print_delta):
        '''
        Auxiliary function to solve the polynomial equation for quasinormal frequencies

        Args:
            d (object): Polynomial in omega / quantization condition
            solver (str): What solver to use (alg / num / mpnum)
            print_delta (bool): Whether to display quantization condition
        '''

        # Display characteristic polynomial if needed
        if (print_delta == True): print(d.expand())

        #   Algebraic equation solver (via sympy)
        if (solver == "alg"): sols = sym.solve(d)
        
        #   Numeric polynomial root solver (via sympy)
        if (solver == "num"): sols = sym.nroots(sym.Poly(d, self.w),n=8, maxsteps=500, cleanup=True)

        # Numerical polynomial root solver (via mpmath, faster and more precise)
        if (solver == "mpnum"):

            # Solver parameters

            # Number of digits to evalf sympy expressions
            dps_evalf = 100
            # Maximum number of steps of solver
            nmax_solve = 10000
            # Extra precision of solver
            xprec_solve = 100
            # Number of digits to print out
            dps_print = 20

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

        # Return solutions to polynomial
        return sols
    

    '''
    AIM ALGORITHM: 
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
        self.l = np.empty(self.n_iter, dtype=object)
        self.s = np.empty(self.n_iter, dtype=object)
        #   Initialize array values
        self.l[0] = self.lambda_0
        self.s[0] = self.s_0
        #   Arrays for lambda' and s' derivative of the parameters
        self.lp = np.empty(self.n_iter, dtype=object)
        self.sp = np.empty(self.n_iter, dtype=object)


    #   Solve via AIM algorithm
    def aim_solve(self, display_all = False, solver = "mpnum", print_delta = False):

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

            # Calculate solutions
            sols = self.aim_poly_solve(d, solver, print_delta)
                    
            #   Display the solution for each iteration
            self.aim_display(which = "AIM", 
                             sols = sols, 
                             display_all = display_all, 
                             n = n)

    '''
    IAIM (improved AIM) ALGORITHM:
    -   Calculate series expansion of lambda_0 and s_0 coefficients
        around the point in which we will be evaluating
    -   Create matrices C, D and store initial coefficients
    -   Compute the necessary coefficients using the recursion relations
    -   From the coefficients in the first row in the two last columns,
        apply the quantization condition
    -   Solve the obtained polynomial in w (omega) and filter solutions
    ''' 


    #   FUNCTION SERIES EXPANSION, REMOVE HIGHER ORDER, RETURN ARRAY WITH COEFFICIENTS
    #   Params: func. a, variable x, around point x0, order N (fixed)
    def iaim_series_coeff(self, a, x, x0):

        coeff = np.empty(self.n_iter, dtype=object)
        for i in range(0, self.n_iter):
            coeff[i] = se.diff(a, x, i).subs(x, x0)/sym.factorial(i)

        return coeff


    def iaim_init(self):

        #   Extra iteration to match AIM's solutions
        self.n_iter = self.n_iter + 1

        #   Matrix of coefficients needed
        self.C = np.zeros((self.n_iter,self.n_iter),dtype=object)
        self.D = np.zeros((self.n_iter,self.n_iter),dtype=object)
        #   First column initialized to coefficients of series expansion
        self.C[:,0] = self.iaim_series_coeff(self.lambda_0,self.x,self.x0)
        self.D[:,0] = self.iaim_series_coeff(self.s_0,self.x,self.x0)
    

    #   Solve via IAIM algorithm
    def iaim_solve(self, solver="mpnum", display_all=False, print_delta=False):

        #   Compute iteratively coefficients / matrix elements
        for n in range(0,self.n_iter-1):
            for i in range(0, self.n_iter - 1 - n):
                self.D[i,n+1] = ( (i+1)*self.D[i+1,n] ).expand()
                self.C[i,n+1] = ( (i+1)*self.C[i+1,n] + self.D[i,n] ).expand()
                for k in range(0, i+1):
                    self.C[i,n+1] = ( self.C[i,n+1] + self.C[k,0]*self.C[i-k,n] ).expand()
                    self.D[i,n+1] = ( self.D[i,n+1] + self.D[k,0]*self.C[i-k,n] ).expand()
   
        #   For each iteration compute the quantization condition and obtain modes
        for n in range(1, self.n_iter-1):

            #   Apply quantization condition and find characteristic polynomial
            d = ( self.D[0,n] * self.C[0,n-1] - self.D[0,n-1] * self.C[0,n] ).expand()

            # Calculate solutions
            sols = self.aim_poly_solve(d, solver, print_delta)

            #   Display the solution for each iteration
            sols_sorted = self.aim_display( which = "IAIM", 
                                            sols = sols, 
                                            display_all = display_all, 
                                            n = n,
                                            n_modes = 15)