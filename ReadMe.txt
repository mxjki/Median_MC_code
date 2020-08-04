#Date: 2020-08-04

Description: this contains the implementations of the proposed method and other compared methods in the simulation study.

########################################
Data (n1=400, n2=400, r=3)
########################################
AU.rds:      		    An n1-by-r matrix, generated from matrix(rnorm(n1*r,0,1),nrow=n1).

AV.rds:      		    An n2-by-r matrix, generated from matrix(rnorm(r*n2,0,1),nrow=n2).

error_cauchy_sim_num_1.rds: An n1-by-n2 error matrix for simulation, generated from matrix(rcauchy(n1 * n2, location = 0, scale = 1),nrow=n1).

omega_sim_num_1.rds:        An n1-by-n2 indicator matrix to denote the missingness of simulation dataset, generated from matrix(rbinom(n1 * n2, 1, 0.2),nrow=n1).

########################################
Auxiliary Funcitons
########################################
SVT:       	  This function performs the singular value soft-thresholding procedures.

KernelK:          This function is the kernel function.

rhotau:           This function calculates the quantile loss. We take $\tau=0.5$ to denote the least absolute deviation loss.

min_abs,update_M: These two functions are auxiliary functions in ADMM algorithm to perform least absolute deviation matrix completion.

LADMCADMM:        This function performs the least absolute deviation matrix completion by ADMM.