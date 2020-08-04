# This is a demo code to illustrate a simulation of
# proposed BLADMC in the paper 'Median Matrix
# Completion: from Embarrassment to Optimality'

# See ReadMe.txt for details of each fitting
# functions.

source("LADMC_functions.R")

####### Data generation

# Generate A_{\star}, here we read the saved data.
AU = readRDS(paste0("AU.rds"))
AV = readRDS(paste0("AV.rds"))
Amat = AU %*% t(AV)

n1 = dim(Amat)[1]  # number of rows of A0
n2 = dim(Amat)[2]  # number of columns of A0

r = dim(AU)[2]  # rank of A0

tau = 0.5  #Denote the median matrix completion


####### Estimation for the proposed BLADMC method.

# Compute BLADMC estimator for simulation 1
errormat = readRDS(paste0("./error_cauchy_sim_num_1.rds"))

omega = readRDS(paste0("./omega_sim_num_1.rds"))

Aobsmat = omega * (Amat + errormat)

# number of row and column blocks
row_block_num = 2
col_block_num = 2

omega_block = omega[1:(n1/row_block_num), 1:(n2/col_block_num)]
mu_block = 1/(2 * sum(omega_block))
maxit_block = 2000
thresh_block = 0.001

trace.it_block = F

lambda_scale

cvLADMC_block = readRDS(paste0("./cvLADMC_block_cauchy_sim_num_1.rds"))

Aini_mat = matrix(0, n1, n2)
for (i in 1:row_block_num) {
    for (j in 1:col_block_num) {
        
        LADMC_block = LADMCADMM(Aobsmat[(1 + (i - 1) * 
            n1/row_block_num):(i * n1/row_block_num), 
            (1 + (j - 1) * n2/col_block_num):(j * n2/col_block_num)], 
            omega[(1 + (i - 1) * n1/row_block_num):(i * 
                n1/row_block_num), (1 + (j - 1) * n2/col_block_num):(j * 
                n2/col_block_num)], n1/row_block_num, 
            n2/col_block_num, tau, mu_block, lambda = svd(Aobsmat[(1 + 
                (i - 1) * n1/row_block_num):(i * n1/row_block_num), 
                (1 + (j - 1) * n2/col_block_num):(j * 
                  n2/col_block_num)])$d[1] * mu_block * 
                cvLADMC_block$lambda_LADMC[cvLADMC_block$cvLADMC_grid], 
            maxit_block, thresh_block, warm.start = NULL, 
            trace.it_block)
        
        Aini_mat[(1 + (i - 1) * n1/row_block_num):(i * 
            n1/row_block_num), (1 + (j - 1) * n2/col_block_num):(j * 
            n2/col_block_num)] = LADMC_block$u %*% (LADMC_block$d * 
            t(LADMC_block$v))
    }
}

saveRDS(Aini_mat, paste0("./BLADMC_cauchy_sim_num_1.rds"))





####### Evaluations for the proposed BLADMC method.

# RMSE_A0
norm(Aini_mat - Amat, type = "F")/sqrt(n1 * n2)

# MAE_A0
mean(abs(Aini_mat - Amat))

# rank
sum(svd(Aini_mat)$d > 1e-10)
