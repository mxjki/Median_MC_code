# This is a demo code to illustrate a simulation of
# proposed DLADMC in the paper 'Median Matrix
# Completion: from Embarrassment to Optimality'

# See ReadMe.txt for details of each fitting
# functions.

source("LADMC_functions.R")
library("softImpute")

####### Data generation

# Generate A_{\star}, here we read the saved data.
AU = readRDS(paste0("AU.rds"))
AV = readRDS(paste0("AV.rds"))
Amat = AU %*% t(AV)

n1 = dim(Amat)[1]  # number of rows of A0
n2 = dim(Amat)[2]  # number of columns of A0

r = dim(AU)[2]  # rank of A0

tau = 0.5  #Denote the median matrix completion

####### Estimation for the proposed DLADMC method.

# Compute DLADMC estimator for simulation 1
errormat = readRDS(paste0("./error_cauchy_sim_num_1.rds"))

omega = readRDS(paste0("./omega_sim_num_1.rds"))

Aobsmat = omega * (Amat + errormat)

Aini_mat = readRDS(paste0("./BLADMC_cauchy_sim_num_1.rds"))

# number of row and column blocks
row_block_num = 2
col_block_num = 2

maxit = 2000
thresh = 1e-05
trace.it = F

N = sum(omega)

hc = 0.1

an0 = hc * sqrt(n1^2 * n2^2 * max(n1/row_block_num, 
    n2/col_block_num) * log(n1/row_block_num + n2/col_block_num)/((n1/row_block_num) * 
    (n2/col_block_num) * N))

for (t in 1:5) {
    
    cvDLADMC = readRDS(paste0("./cvDLADMC_cauchy_iter_", 
        t, ".rds"))
    
    h = hc * (sqrt(r * n1 * n2 * max(n1, n2) * log(n1 + 
        n2)/N) + r^(-1/2) * min(n1, n2) * (an0 * r^(1/2)/min(n1, 
        n2))^(2^t))/sqrt(n1 * n2)
    
    Kernelmat = apply((Aobsmat - Aini_mat)/h, c(1, 2), 
        KernelK)
    
    fhat = sum(omega * Kernelmat/h)/N
    
    Aobsmat_tilde = Aini_mat - (matrix(as.numeric(Aobsmat - 
        Aini_mat <= 0), n1, n2) - tau)/fhat
    
    AobsDLADMC = omega * Aobsmat_tilde
    
    AobsDLADMC_NA = AobsDLADMC
    AobsDLADMC_NA[AobsDLADMC_NA == 0] = NA
    svdAobsDLADMC_NA = svd(AobsDLADMC)
    
    DLADMC = softImpute(AobsDLADMC_NA, rank.max = min(n1, 
        n2) - 1, lambda = svdAobsDLADMC_NA$d[1] * cvDLADMC$lambda_MHT[cvDLADMC$cvMHT_grid], 
        type = c("svd"), thresh, warm.start = NULL, 
        maxit, trace.it)
    
    Aini_mat = DLADMC$u %*% (DLADMC$d * t(DLADMC$v))
    
}

saveRDS(DLADMC, paste0("./DLADMC_cauchy_sim_num_1.rds"))



####### Evaluations for the proposed DLADMC method.

# RMSE_A0
norm(DLADMC$u %*% (DLADMC$d * t(DLADMC$v)) - Amat, type = "F")/sqrt(n1 * 
    n2)

# MAE_A0
mean(abs(DLADMC$u %*% (DLADMC$d * t(DLADMC$v)) - Amat))

# rank
sum(DLADMC$d > 1e-10)
