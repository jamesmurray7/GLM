rm(list=ls())
# Determine OS and setup directory to save in.
save.dir <- unname(ifelse(Sys.info()[1]=='Linux', '/data/c0061461/cmp-grids/', './data.nosync/'))
sourceCpp('grid-test.cpp')

# \lambda grid ----
lambda.mat <- gen_lambda_mat(1000,10)
# Compare this to Pete's grid
load(paste0(save.dir, 'grid_1K.RData'))

dim(lambda.mat);dim(grid_1K)
lambda.mat[999,]; grid_1K[999,]
which(lambda.mat[999,] == 1e3)
which(grid_1K[999,] == 0)

sapply(1:999, function(i) min(abs(grid_1K[i,]-lambda.mat[i,]))) # In all cases, there is at least some agreement
sapply(1:999, function(i) max(abs(grid_1K[i,]-lambda.mat[i,]))) # Appears to start to diverge at 152nd row of mu.

# Grids for variances and logZ ----
V.mat <- gen_V_mat(1000, 100, lambda.mat)
logZ.mat <- gen_logZ_mat(1000, 100, lambda.mat)

# Save all ----
save(lambda.mat, file = paste0(save.dir, 'lambda.RData'))
save(V.mat, file = paste0(save.dir, 'V.RData'))
save(logZ.mat, file = paste0(save.dir, 'logZ.RData'))

# Generate these using Pete's lambda grid instead? ----
V.mat <- gen_V_mat(1000, 100, grid_1K)
logZ.mat <- gen_logZ_mat(1000, 100, grid_1K)
save(V.mat, file =  paste0(save.dir, 'V_Pete.RData'))
save(logZ.mat, file =  paste0(save.dir, 'logZ_Pete.RData'))


# 10K version -------------------------------------------------------------
lambda.mat <- gen_lambda_mat(10000, 10)
save(lambda.mat, file = paste0(save.dir, 'lambda10K.RData'))
logZ.mat <- gen_logZ_mat(10000, 100, lambda.mat)
save(logZ.mat, file = paste0(save.dir, 'logZ10K.RData'))
V.mat <- gen_V_mat(10000, 100, lambda.mat)
save(V.mat, file = paste0(save.dir, 'V10K.RData'))

# Generate these using Pete's lambda grid instead? ----
load(paste0(save.dir, 'lambda10K_Pete.RData'))
V.mat <- gen_V_mat(10000, 100, lam_grid_poly_10K[-10000,])
save(V.mat, file =  paste0(save.dir, 'V10K_Pete.RData'))
logZ.mat <- gen_logZ_mat(10000, 100, lam_grid_poly_10K[-10000,])
save(logZ.mat, file =  paste0(save.dir, 'logZ10K_Pete.RData'))
