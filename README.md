# TFG-QNMs-BH
Numerical code used for the Undergraduate Dissertation: Quasinormal Modes in Black Holes

Implementation of Asymptotic Iteration Method found in main code "aim_qnm.py", which is called by the script "aim_qnm_(...).y" found for each case in the subfolders "schwarzschild" and "poschl-teller". We must specify as inputs functions lambda_0 and s_0 previously computed (evaluating all parameters numerically except the coordinate and eigenavlue), the evaluation point, and the number of iterations to perform. 

The jupyter notebook files include a brief description on the calculations done for the Poschl-Teller case, which explains the idea behind the AIM method.

All used scripts to generate the graphs can be found in the folder "plots", while the data obtained while looping over position of momentum can be found in subfolder "loops". 
