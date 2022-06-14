rm(list=ls())
# Determine OS and setup directory to save in.
save.dir <- unname(ifelse(Sys.info()[1]=='Linux', '/data/c0061461/cmp-grids/', './data.nosync/'))
sourceCpp('grid-test.cpp')

# \lambda grid ----
lambda.mat <- gen_lambda_mat(1000,10)
# Compare this to Pete's grid
load(paste0(save.dir, 'grid_1K.RData'))

comp <- lambda.mat[-nrow(lambda.mat),]
dim(comp);dim(grid_1K)
comp[999,]; grid_1K[999,]

sapply(1:999, function(i) min(abs(grid_1K[i,]-comp[i,]))) # In all cases, there is at least some agreement
sapply(1:999, function(i) max(abs(grid_1K[i,]-comp[i,]))) # Appears to start to diverge at 152nd row of mu.

# Grids for variances and logZ ----
V.mat <- gen_V_mat(1000, 100, lambda.mat)
logZ.mat <- gen_logZ_mat(1000, 100, lambda.mat)


# Save all ----
save(lambda.mat, file = paste0(save.dir, 'lambda.RData'))
save(V.mat, file = paste0(save.dir, 'V.RData'))
save(logZ.mat, file = paste0(save.dir, 'logZ.RData'))
