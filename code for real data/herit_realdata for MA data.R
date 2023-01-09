
library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)
library(RSpectra)
noise = 1

load('/data2/thetm/MA/new.masnps.rda')
masnps.filter = masnps.filter/2

p = ncol(masnps.filter)
n = nrow(masnps.filter)
source('/data2/thetm/review_HERI/mainfunctions.R')


# convex optimization method
start_time <- Sys.time()
eprsm <- EigenPrism(ypenicilin, masnps.filter, target = 'heritability')
end_time <- Sys.time()
end_time - start_time

# MLE method
start_time <- Sys.time()
M = tcrossprod(masnps.filter)
O = eigen(M)$vectors
lambda = eigen(M)$values
yTilde = t(O) %*% ypenicilin
mle.r =ridge.mle(yTilde,lambda,1.e-6)
eta_chap = 1 - mle.r$sig2hat / var(ypenicilin)
s_eta = mle.r$sig2hat/sqrt(2)
eta_chap + qnorm(1-0.05/2)* s_eta*c(-1,1) 
end_time <- Sys.time()
end_time - start_time

# moment method
start_time <- Sys.time()
mmdker = mm2dicker(ypenicilin, masnps.filter,CI = T )
end_time <- Sys.time()
end_time - start_time

# scaled lasso method
start_time <- Sys.time()
scaled_lasso <- scaled.lasso(masnps.filter , ypenicilin)
Slasso = 1-scaled_lasso$sigma2^2/var(ypenicilin)
k = scaled_lasso$non.zeros
#tam = svd(masnps.filter,0,0)
Slasso- c(-1,1)*log(.5)*(k*log(p)/n + 1/sqrt(n))
end_time <- Sys.time()
end_time - start_time

# Enet method
start_time <- Sys.time()
glmnetcv = cv.glmnet(masnps.filter, ypenicilin, alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
enet = coef( glmnetcv,  s='lambda.min')[-1]
Enet = enet[tam] %*% cov2(masnps.filter[,tam]) %*% enet[tam]/ var(ypenicilin)
end_time <- Sys.time()
end_time - start_time
