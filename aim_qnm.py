#   Importing required libraries
import numpy as np      #   Used for numerical calculations
import sympy as sym     #   Used for symbolic calculations: Series expansion,
import time             #   Used to calculate computation times in critical parts
                        #   (better to do it in main code)

#   Supposedly faster/more precise, unused for now
import symengine as se
import flint as fl

'''
Improved AIM algorithm (IAIM) taylored for computing quasinormal modes (QNMs)
of black holes (BHs).
Focused on solving second-order Schrödinger-like equation already written in the "AIM" form:

d^2/dx^2 psi(x) = lambda_0(x) d/dx psi(x) + s_0(x) psi(x)

Therefore, we have to compute the parameters lambda_0 and s_0 beforehand manually
(rewriting the potential in appropriate coordinates and including explicitly boundary conditions
into the solution, that is, incoming waves at horizon and outgoing at infinity)

(Possible improvement: code to compute lambda_0 and s_0 for any kind of potential)

Inputs / Initial parameters:

    -   lambda_0, s_0 - > Expressions of initial parameters lambda_0 and s_0 in terms of only numeric values and symbolic coordinate (x)
    -   x_0           - > Point of evaluation (x0) around which to perform series expansion / substitutions

    ~   (Perform all substitutions in main code) Numeric values of all parameters included in lambda_0 and s_0 (such as mass, spin, angular momentum, etc.)

Options:

    -   re_sign, im_sign - > What solutions to filter (positive real part, negative imaginary, etc.)
    -   alg = sym/num    - > What algorithm to use to solve the polynomial (algebraic/symbolic or numeric solutions)  
    -   n_iter           - > Number of iterations to perform
    -   prec             - > Decimal places to display

Outputs:

    -   Array of all filtered/unfiltered solutions (complex frequencies omega) ordered by mode number (n)
    -   Total computation time
    -   Text update that ensures program is not stuck in loop
    -   Optional additional parameters/values of interest
    -   Write to file

Algorithms implemented:

    * Asymptotic Iteration Method (AIM)
    
    * Improves Asymptotic Iteration Method (IAIM)

Tested for:

    -   ...

'''

class aim_solver(object):

    w = sym.symbols("\omega") #  Symbol representing complex frequency of qnms

    #   Create AIM object with the necessary initial parameters, and create necessary arrays for lambda_n, s_n and derivatives
    #   Initial parameters depend only on symbols (x,y,...) defined in main code as necessary, with all numeric substitutions
    #   already performed

    def __init__(self, lambda_0, s_0, n_iter): 

        
        #   Define initial parameters
        self.lambda_0 = lambda_0
        self.s_0 = s_0
        #   Number of iterations (size of arrays / matrices) to perform
        self.n_iter = n_iter

        return
    
    '''
    AIM algorithm: 
    -   Begin by defining initial parameters lambda_0 and s_0
    -   Calculate necessary derivatives of "previous" parameters
        lambda'_(n-1) and s'_(n-1)
    -   Calculate new parameters lambda_n and s_n in terms of previous
        parameters and their derivatives
    -   Obtain delta (quantization condition)
    -   Evaluate at point of interest
    -   Solve polynomial in w (omega) algebraically and filter solutions
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

        return
    
    #   Function to display results from AIM algorithm
    def aim_display(self, sols, display_all):

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
            print("w_" + str(i) + " = " +  str(sols_sorted[i]))

        return
    
    #   Function that performs the AIM algorithm, calculating all parameters lambda_n, s_n and derivatives
    #   Params: display (True/False) -> Shows solution for each step
    #           x, x0 (symbol, num)  -> Variable to differentiate and evaluation point

    def aim_alg_solve(self, x, x0, display_all):

        for n in range(1, self.n_iter):

            print("\n*** ITERATION n=" + str(n) + " ***\n")

            #   Calculate the previous derivatives
            self.lp[n-1] = sym.diff(self.l[n-1],x)
            self.sp[n-1] = sym.diff(self.s[n-1],x)

            #   Calculate the new parameters lambda_n and s_n
            #   Using the previous lambda_n-1, s_n-1, and derivatives lambda'_n-1 and s'_n-1
            self.l[n] = self.lp[n-1] + self.s[n-1] + self.l[0]*self.l[n-1]
            self.s[n] = self.sp[n-1] + self.s[0]*self.l[n-1]

            #   Quantization condition delta
            d = self.s[n]*self.l[n-1] - self.s[n-1]*self.l[n]

            #   "Characteristic polynomial" obtained by evaluating quantization condition
            p = d.subs(x,x0)

            #   Algebraic equation solver (via sympy)
            sols = sym.solve(p) #   Solve the characteristic polynomial
            self.aim_display(sols, display_all) #   Display the solution for each iteration
            
    def aim_num_solve(self, x, x0, display_all):

        for n in range(1, self.n_iter):

            print("\n*** ITERATION n=" + str(n) + " ***\n")

            #   Calculate the previous derivatives
            self.lp[n-1] = sym.diff(self.l[n-1],x)
            self.sp[n-1] = sym.diff(self.s[n-1],x)

            #   Calculate the new parameters lambda_n and s_n
            #   Using the previous lambda_n-1, s_n-1, and derivatives lambda'_n-1 and s'_n-1
            self.l[n] = self.lp[n-1] + self.s[n-1] + self.l[0]*self.l[n-1]
            self.s[n] = self.sp[n-1] + self.s[0]*self.l[n-1]

            #   Quantization condition delta
            d = self.s[n]*self.l[n-1] - self.s[n-1]*self.l[n]

            #   "Characteristic polynomial" obtained by evaluating quantization condition
            p = d.subs(x,x0)

            #   Construct a sympy polynomial
            ppol = sym.Poly(p,self.w)
            
            #   Numeric polynomial root solver (via sympy)
            sols = sym.nroots(ppol) #   Find roots of characteristic polynomial
            self.aim_display(sols, display_all) #   Display the solution for each iteration

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