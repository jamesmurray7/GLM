rm(list=ls())
# Determine OS and setup directory to save in.
save.dir <- unname(ifelse(Sys.info()[1]=='Linux', '/data/c0061461/cmp-grids/new', paste0(getwd(),'/data.nosync/')))
sourceCpp('testing.cpp')

# 10K version -------------------------------------------------------------
lambda.mat <- gen_lambda_mat(10000, 10)
save(lambda.mat, file = paste0(save.dir, 'lambda10K.RData'))
logZ.mat <- gen_logZ_mat(10000, 100, lambda.mat)
save(logZ.mat, file = paste0(save.dir, 'logZ10K.RData'))
V.mat <- gen_V_mat(10000, 100, lambda.mat, logZ.mat)
save(V.mat, file = paste0(save.dir, 'V10K.RData'))

# # Generate these using Pete's lambda grid instead? ----
# load(paste0(save.dir, 'lambda10K_Pete.RData'))
# logZ.mat <- gen_logZ_mat(10000, 100, lam_grid_poly_10K[-10000,])
# save(logZ.mat, file =  paste0(save.dir, 'logZ10K_Pete.RData'))
# V.mat <- gen_V_mat(10000, 100, lam_grid_poly_10K[-10000,], logZ.mat)
# save(V.mat, file =  paste0(save.dir, 'V10K_Pete.RData'))
