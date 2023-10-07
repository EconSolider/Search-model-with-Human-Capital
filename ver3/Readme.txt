Program of search model with human capital accumulation and 
	variable search intensity

Main Programs:
runner.m  : main file to calculate steady state and run dynare files
a_param_model.mod : dyanre file to contain the main model

b_analysisA_stoc: dyanre file to perform stochastic simulation of exogenous shocks

c_analysisB_deter: dynare file to perform determinstic simulation
(model's reaction to permanant change in search and learning time, and search time)

Functions:
main_fun_a: function to calibrate z, phi, lambda
main_fun_b: function to solve value W and U of workers

Insert Programs (trivial)
Set_param_value.inc : dynare insert file to set parameter value
b_analysisA_stoc_steadystate.m: solve steady state for analysis A.
c_analysisB_deter_steadystate.m: solve steady state for analysis B.
